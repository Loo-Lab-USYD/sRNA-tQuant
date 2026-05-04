#!/usr/bin/env nextflow

/*
 * tRNA Isodecoder Quantification Pipeline
 *
 * Quantify cell-line and tissue-specific tRNA abundances at anticodon-level
 * resolution from small RNA-seq data.
 *
 * Two quantification approaches (configurable via params.quantifier):
 *   - 'salmon'   (default): probabilistic EM via Salmon
 *   - 'decision' (legacy):  anticodon-aware binary keep/discard
 */

nextflow.enable.dsl = 2

// Module imports
include { TRNASCAN_SE; MASK_GENOME; BUILD_PRETRNA; BUILD_ARTIFICIAL_GENOME;
          BUILD_MATURE_LIBRARY; INDEX_REFERENCES } from './modules/prepare_reference'
include { TRIM_READS }                            from './modules/trim_reads'
include { REMOVE_RRNA }                           from './modules/remove_rrna'
include { MAP_ARTIFICIAL }                        from './modules/map_artificial'
include { MAP_CLUSTERS }                          from './modules/map_clusters'
include { FILTER_TOPSCORE; SALMON_QUANT;
          SALMON_TO_RPM }                         from './modules/quantify_salmon'
include { RESOLVE_MULTIMAPPERS; QUANTIFY_RPM }    from './modules/quantify_decision'
include { CALL_VARIANTS; COUNT_ALLELES;
          ADJUST_MODIFICATIONS }                  from './modules/call_modifications'
include { AGGREGATE_RESULTS }            from './modules/aggregate'

/*
 * Validate parameters
 */
def validateParams() {
    if (!params.input) {
        error "Please provide a samplesheet via --input"
    }
    if (!params.genome_fasta && !params.trna_reference) {
        error "Please provide either --genome_fasta or --trna_reference"
    }
    if (params.genome_fasta && params.trna_reference) {
        error "Please provide only one of --genome_fasta or --trna_reference, not both"
    }
    if (!(params.quantifier in ['salmon', 'decision'])) {
        error "params.quantifier must be 'salmon' or 'decision', got: ${params.quantifier}"
    }
}

/*
 * Parse samplesheet CSV into channel of [meta, reads] tuples
 * Expected columns: sample_id, fastq_1, fastq_2 (optional)
 */
def parseSamplesheet(samplesheet) {
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id, single_end: !row.fastq_2]
            def reads = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)]
                                    : file(row.fastq_1)
            [meta, reads]
        }
}

/*
 * Sub-workflow: prepare reference from genome FASTA
 */
workflow PREPARE_REFERENCE {
    take:
    genome_fasta

    main:
    TRNASCAN_SE(genome_fasta)
    MASK_GENOME(genome_fasta, TRNASCAN_SE.out.bed12)
    BUILD_PRETRNA(genome_fasta, TRNASCAN_SE.out.bed12)
    BUILD_ARTIFICIAL_GENOME(MASK_GENOME.out.masked_fasta, BUILD_PRETRNA.out.pretrna_fasta)
    BUILD_MATURE_LIBRARY(genome_fasta, TRNASCAN_SE.out.bed12)
    INDEX_REFERENCES(BUILD_ARTIFICIAL_GENOME.out.fasta, BUILD_MATURE_LIBRARY.out.cluster_fasta)

    emit:
    artificial_index = INDEX_REFERENCES.out.artificial_index
    cluster_index    = INDEX_REFERENCES.out.cluster_index
    cluster_fasta    = BUILD_MATURE_LIBRARY.out.cluster_fasta
    cluster_info     = BUILD_MATURE_LIBRARY.out.cluster_info
    cluster_fai      = INDEX_REFERENCES.out.cluster_fai
    trna_list        = BUILD_MATURE_LIBRARY.out.trna_list
    pretrna_fasta    = BUILD_PRETRNA.out.pretrna_fasta
    pretrna_bed12    = BUILD_PRETRNA.out.pretrna_bed12
}

/*
 * Main workflow
 */
workflow {
    validateParams()

    // Parse input samples
    ch_samples = parseSamplesheet(params.input)

    // --- Reference preparation ---
    if (params.genome_fasta) {
        PREPARE_REFERENCE(file(params.genome_fasta))
        ch_artificial_index = PREPARE_REFERENCE.out.artificial_index
        ch_cluster_index    = PREPARE_REFERENCE.out.cluster_index
        ch_cluster_fasta    = PREPARE_REFERENCE.out.cluster_fasta
        ch_cluster_info     = PREPARE_REFERENCE.out.cluster_info
        ch_cluster_fai      = PREPARE_REFERENCE.out.cluster_fai
        ch_trna_list        = PREPARE_REFERENCE.out.trna_list
        ch_pretrna_fasta    = PREPARE_REFERENCE.out.pretrna_fasta
        ch_pretrna_bed12    = PREPARE_REFERENCE.out.pretrna_bed12
    } else {
        // Pre-computed reference bundle
        def ref = params.trna_reference
        ch_artificial_index = file("${ref}/artificial_genome_index")
        ch_cluster_index    = file("${ref}/cluster_index")
        ch_cluster_fasta    = file("${ref}/cluster.fa")
        ch_cluster_info     = file("${ref}/clusterInfo.fa")
        ch_cluster_fai      = file("${ref}/cluster.fa.fai")
        ch_trna_list        = file("${ref}/tRNAs.txt")
        ch_pretrna_fasta    = file("${ref}/pre_trnas.fa")
        ch_pretrna_bed12    = file("${ref}/pre_trnas.bed12")
    }

    // --- Per-sample processing ---

    // 1. Trim reads
    TRIM_READS(ch_samples)

    // 1b. Remove rRNA (optional — provide --rrna_index to enable)
    ch_trimmed = TRIM_READS.out.reads
    if (params.rrna_index) {
        REMOVE_RRNA(ch_trimmed, file(params.rrna_index))
        ch_post_rrna = REMOVE_RRNA.out.reads
    } else {
        ch_post_rrna = ch_trimmed
    }

    // 2-4. Map to artificial genome, filter genome reads, filter mature reads
    MAP_ARTIFICIAL(
        ch_post_rrna,
        ch_artificial_index,
        ch_pretrna_fasta,
        ch_pretrna_bed12
    )

    // 5. Map to clustered tRNA reference
    MAP_CLUSTERS(
        MAP_ARTIFICIAL.out.reads,
        ch_cluster_index,
        params.quantifier
    )

    // 6-7. Quantify (Salmon or Decision path)
    if (params.quantifier == 'salmon') {
        FILTER_TOPSCORE(MAP_CLUSTERS.out.sam)
        SALMON_QUANT(FILTER_TOPSCORE.out.bam, ch_cluster_fasta)

        // Join Salmon output with topscore BAM by sample ID for span-weighted mode
        ch_quant_bam = SALMON_QUANT.out.quant_dir.join(FILTER_TOPSCORE.out.bam)
        SALMON_TO_RPM(
            ch_quant_bam.map { meta, qdir, bam -> [meta, qdir] },
            ch_cluster_info,
            ch_quant_bam.map { meta, qdir, bam -> [meta, bam] }
        )
        ch_rpm = SALMON_TO_RPM.out.rpm
        ch_cluster_rpm = SALMON_TO_RPM.out.cluster_rpm
    } else {
        RESOLVE_MULTIMAPPERS(MAP_CLUSTERS.out.sam, ch_cluster_info)
        QUANTIFY_RPM(RESOLVE_MULTIMAPPERS.out.sam, ch_cluster_info, ch_trna_list)
        ch_rpm = QUANTIFY_RPM.out.rpm
        ch_cluster_rpm = QUANTIFY_RPM.out.cluster_rpm
    }

    // 8-9. Modification calling (optional)
    if (!params.skip_modification_calling) {
        // Use unique-only reads for variant calling
        CALL_VARIANTS(MAP_CLUSTERS.out.sam, ch_cluster_fasta, ch_cluster_fai)

        // Join SAM and VCF by sample ID to prevent sample mismatches
        ch_sam_vcf = MAP_CLUSTERS.out.sam.join(CALL_VARIANTS.out.vcf)
        COUNT_ALLELES(
            ch_sam_vcf.map { meta, sam, vcf -> [meta, sam] },
            ch_sam_vcf.map { meta, sam, vcf -> [meta, vcf] },
            ch_cluster_fasta, ch_cluster_fai
        )

        // Join RPM, ASE, and cluster RPM by sample ID
        ch_rpm_ase_cr = ch_rpm.join(COUNT_ALLELES.out.ase).join(ch_cluster_rpm)
        ADJUST_MODIFICATIONS(
            ch_rpm_ase_cr.map { meta, rpm, ase, cr -> [meta, rpm] },
            ch_rpm_ase_cr.map { meta, rpm, ase, cr -> [meta, ase] },
            ch_cluster_info, ch_trna_list,
            ch_rpm_ase_cr.map { meta, rpm, ase, cr -> [meta, cr] }
        )
        ch_final_rpm = ADJUST_MODIFICATIONS.out.rpm
        ch_ase_files = COUNT_ALLELES.out.ase.map { meta, ase -> ase }.collect()
    } else {
        ch_final_rpm = ch_rpm
        ch_ase_files = file('NO_ASE')
    }

    // --- Aggregation ---
    ch_all_rpms = ch_final_rpm.map { meta, rpm -> rpm }.collect()

    AGGREGATE_RESULTS(ch_all_rpms, ch_ase_files)

    // MultiQC
    // MULTIQC disabled: ch_multiqc = Channel.empty()
    // MULTIQC disabled: .mix(TRIM_READS.out.json.collect())
    // MULTIQC disabled: .mix(MAP_CLUSTERS.out.log.collect())
    // MULTIQC disabled
}

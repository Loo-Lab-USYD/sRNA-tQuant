/*
 * PREPARE_REFERENCE: predict tRNA genes, mask genome, build artificial genome,
 * create mature tRNA library, cluster, and index.
 *
 * Only runs when params.genome_fasta is provided (no pre-computed bundle).
 */

process TRNASCAN_SE {
    tag "tRNAscan-SE"
    label 'process_high'

    input:
    path genome_fasta

    output:
    path "trnascan_nuclear.bed12", emit: nuclear_bed
    path "trnascan_mito.bed12",    emit: mito_bed
    path "trnascan_merged.bed12",  emit: bed12

    script:
    """
    # Nuclear tRNAs
    tRNAscan-SE -q --thread ${task.cpus} -b trnascan_nuclear.bed12 ${genome_fasta}

    # Mitochondrial tRNAs (organellar mode)
    samtools faidx ${genome_fasta}
    # Extract chrM/MT sequence
    if samtools faidx ${genome_fasta} chrM > chrM.fa 2>/dev/null; then
        tRNAscan-SE -q -O -b trnascan_mito.bed12 chrM.fa
    elif samtools faidx ${genome_fasta} MT > chrM.fa 2>/dev/null; then
        tRNAscan-SE -q -O -b trnascan_mito.bed12 chrM.fa
    else
        touch trnascan_mito.bed12
    fi

    # Merge (exclude chrM from nuclear to avoid duplicates)
    grep -v -E '^chrM|^MT' trnascan_nuclear.bed12 > trnascan_nuclear_noM.bed12 || true
    cat trnascan_nuclear_noM.bed12 trnascan_mito.bed12 > trnascan_merged.bed12
    """
}

process MASK_GENOME {
    tag "mask genome"
    label 'process_medium'

    input:
    path genome_fasta
    path trna_bed12

    output:
    path "genome_masked.fa", emit: masked_fasta

    script:
    """
    bedtools maskfasta -fi ${genome_fasta} -fo genome_masked.fa -mc N -bed ${trna_bed12}
    """
}

process BUILD_PRETRNA {
    tag "build pre-tRNAs"
    label 'process_low'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path genome_fasta
    path trna_bed12

    output:
    path "pre_trnas.bed12", emit: pretrna_bed12
    path "pre_trnas.fa",    emit: pretrna_fasta

    script:
    """
    # Add 50 nt flanking regions (v2: fix block sizes for -split extraction)
    echo "mod_bed12 v2" > /dev/null
    python3 ${projectDir}/bin/mod_bed12.py ${trna_bed12} pre_trnas.bed12

    # Extract sequences (introns removed via -split)
    bedtools getfasta -name -split -s -fi ${genome_fasta} -bed pre_trnas.bed12 -fo pre_trnas.fa
    """
}

process BUILD_ARTIFICIAL_GENOME {
    tag "artificial genome"
    label 'process_low'

    input:
    path masked_fasta
    path pretrna_fasta

    output:
    path "artificial_genome.fa", emit: fasta

    script:
    """
    cat ${masked_fasta} ${pretrna_fasta} > artificial_genome.fa
    """
}

process BUILD_MATURE_LIBRARY {
    tag "mature tRNAs"
    label 'process_low'
    publishDir "${params.outdir}/reference", mode: 'copy',
        pattern: '{cluster.fa,clusterInfo.fa,tRNAs.txt}'

    input:
    path genome_fasta
    path trna_bed12

    output:
    path "mature_trnas.fa",    emit: mature_fasta
    path "cluster.fa",         emit: cluster_fasta
    path "clusterInfo.fa",     emit: cluster_info
    path "tRNAs.txt",          emit: trna_list

    script:
    """
    # Extract mature tRNA sequences
    bedtools getfasta -name -split -s -fi ${genome_fasta} -bed ${trna_bed12} -fo raw_trnas.fa

    # Add CCA tail and filter pseudogenes
    add_cca.py raw_trnas.fa mature_trnas.fa

    # Cluster identical sequences
    cluster_trnas.py mature_trnas.fa cluster.fa clusterInfo.fa

    # Generate tRNA anticodon list from clusterInfo
    grep '^>' clusterInfo.fa | sed 's/.*-\\(i\\?[A-Z][a-zC]*[ACTGN]\\{3\\}\\)[:(].*/\\1/' | sort -u > tRNAs.txt
    """
}

process INDEX_REFERENCES {
    tag "index"
    label 'process_medium'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path artificial_fasta
    path cluster_fasta

    output:
    path "artificial_genome_index", emit: artificial_index
    path "cluster_index",           emit: cluster_index
    path "${cluster_fasta}.fai",    emit: cluster_fai

    script:
    """
    # Index artificial genome for bowtie2
    mkdir artificial_genome_index
    bowtie2-build ${artificial_fasta} artificial_genome_index/artificial_genome

    # Index clustered tRNAs for bowtie2
    mkdir cluster_index
    bowtie2-build ${cluster_fasta} cluster_index/cluster

    # samtools index for variant calling
    samtools faidx ${cluster_fasta}
    """
}

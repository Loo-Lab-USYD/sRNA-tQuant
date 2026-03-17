/*
 * QUANTIFY_SALMON: filter top-scoring alignments and run Salmon in
 * alignment-based mode for probabilistic multimapper resolution.
 */

process FILTER_TOPSCORE {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}_topscore.bam"), emit: bam

    script:
    """
    # Filter to equal top-scoring alignments
    filter_topscore.py ${sam} ${meta.id}_topscore.sam

    echo "filter_topscore v2" > /dev/null

    # Convert to BAM (Salmon accepts BAM; must NOT be position-sorted)
    samtools view -@ ${task.cpus} -bS ${meta.id}_topscore.sam > ${meta.id}_topscore.bam
    rm -f ${meta.id}_topscore.sam
    """
}

process SALMON_QUANT {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(bam)
    path cluster_fasta

    output:
    tuple val(meta), path("${meta.id}_salmon"), emit: quant_dir
    path "${meta.id}_salmon/quant.sf",          emit: quant_sf

    script:
    def extra_args = task.ext.args ?: ''
    """
    salmon quant \\
        --libType A \\
        --alignments ${bam} \\
        --targets ${cluster_fasta} \\
        --output ${meta.id}_salmon \\
        --threads ${task.cpus} \\
        ${extra_args}
    """
}

process SALMON_TO_RPM {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(quant_dir)
    path cluster_info
    tuple val(meta2), path(bam)

    output:
    tuple val(meta), path("${meta.id}.RPM.csv"),         emit: rpm
    tuple val(meta), path("${meta.id}.cluster_rpm.csv"), emit: cluster_rpm
    tuple val(meta), path("${meta.id}.RPM_annotations.csv"), emit: annotations, optional: true

    script:
    def bam_arg = params.span_weighted ? "--bam ${bam} --annotations ${meta.id}.RPM_annotations.csv" : ''
    """
    echo "salmon_to_rpm v3" > /dev/null
    span_weighted_rpm.py ${quant_dir}/quant.sf ${cluster_info} ${meta.id}.RPM.csv \\
        --cluster-rpm ${meta.id}.cluster_rpm.csv \\
        ${bam_arg}
    """
}

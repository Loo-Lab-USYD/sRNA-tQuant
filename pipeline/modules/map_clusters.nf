/*
 * MAP_CLUSTERS: map mature tRNA reads to the clustered tRNA reference.
 * Bowtie2 flags depend on which quantifier is selected downstream.
 */

process MAP_CLUSTERS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    path cluster_index
    val quantifier

    output:
    tuple val(meta), path("${meta.id}_clusters.sam"), emit: sam
    path "${meta.id}_cluster_mapping.log",            emit: log

    script:
    // Salmon needs -a (all alignments); Decision needs -k 100
    def multimapper_flag = quantifier == 'salmon' ? '-a' : '-k 100'
    def extra_args = task.ext.args ?: ''
    """
    bowtie2 \\
        --very-sensitive \\
        -p ${task.cpus} \\
        -L 14 -D 25 \\
        ${multimapper_flag} \\
        --no-unal \\
        -x ${cluster_index}/cluster \\
        -U ${reads} \\
        -S ${meta.id}_clusters.sam \\
        ${extra_args} \\
        2> ${meta.id}_cluster_mapping.log
    """
}

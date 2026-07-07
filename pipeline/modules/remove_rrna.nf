/*
 * REMOVE_RRNA: deplete rRNA reads using bowtie2 --un-gz.
 * Reads that do NOT map to the rRNA index are passed downstream.
 * Only runs if --rrna_index is provided.
 */

process REMOVE_RRNA {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    path rrna_index

    output:
    tuple val(meta), path("${meta.id}_no_rrna.fastq.gz"), emit: reads
    path "${meta.id}_rrna.log",                            emit: log

    script:
    """
    bowtie2 \\
        -x ${rrna_index}/rrna \\
        -U ${reads} \\
        --un-gz ${meta.id}_no_rrna.fastq.gz \\
        --very-fast \\
        -p ${task.cpus} \\
        -S /dev/null \\
        2> ${meta.id}_rrna.log
    """
}

/*
 * TRIM_READS: adapter and quality trimming with fastp.
 * Replaces BBduk. Built-in QC report (MultiQC-compatible).
 */

process TRIM_READS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: reads
    path "${meta.id}_fastp.json",                         emit: json
    path "${meta.id}_fastp.html",                         emit: html

    script:
    def args = task.ext.args ?: ''
    def adapter_arg = params.adapter_sequence ? "--adapter_sequence ${params.adapter_sequence}" : ''
    if (meta.single_end) {
        """
        fastp \\
            -i ${reads} \\
            -o ${meta.id}_trimmed.fastq.gz \\
            --json ${meta.id}_fastp.json \\
            --html ${meta.id}_fastp.html \\
            --qualified_quality_phred 25 \\
            --length_required ${params.min_length} \\
            --length_limit ${params.max_length_se} \\
            ${adapter_arg} \\
            --thread ${task.cpus} \\
            ${args}
        """
    } else {
        """
        fastp \\
            -i ${reads[0]} \\
            -I ${reads[1]} \\
            -o ${meta.id}_trimmed_R1.fastq.gz \\
            -O ${meta.id}_trimmed_R2.fastq.gz \\
            --json ${meta.id}_fastp.json \\
            --html ${meta.id}_fastp.html \\
            --qualified_quality_phred 25 \\
            --length_required ${params.min_length} \\
            --length_limit ${params.max_length_pe} \\
            ${adapter_arg} \\
            --thread ${task.cpus} \\
            ${args}
        """
    }
}

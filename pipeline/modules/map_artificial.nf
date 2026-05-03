/*
 * MAP_ARTIFICIAL: map trimmed reads to the artificial genome (masked genome +
 * pre-tRNAs) using Bowtie2, then filter to isolate mature tRNA reads.
 */

process MAP_ARTIFICIAL {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    path artificial_index
    path pretrna_fasta
    path pretrna_bed12

    output:
    tuple val(meta), path("${meta.id}_mature.fastq.gz"), emit: reads
    path "${meta.id}_mapping.log",                       emit: log

    script:
    """
    # Map to artificial genome
    bowtie2 \\
        --very-sensitive \\
        -p ${task.cpus} \\
        -L 14 -D 25 \\
        -k 100 \\
        --no-unal \\
        -x ${artificial_index}/artificial_genome \\
        -U ${reads} \\
        -S ${meta.id}_artificial.sam \\
        2> ${meta.id}_mapping.log

    # Filter: remove reads that map to genome (non-pre-tRNA contigs)
    echo "filter_genome_reads v3" > /dev/null
    filter_genome_reads.py \\
        ${pretrna_fasta} \\
        ${meta.id}_artificial.sam \\
        ${meta.id}_filtered.sam \\
        --threads ${task.cpus}

    # Filter: keep only reads mapping to mature tRNA region
    echo "filter_mature_reads v4" > /dev/null
    filter_mature_reads.py \\
        ${pretrna_bed12} \\
        ${meta.id}_filtered.sam \\
        ${reads} \\
        --threads ${task.cpus} \\
        > ${meta.id}_mature.fastq

    # Check for low tRNA mapping rate
    TRNA_READS=\$(( \$(wc -l < ${meta.id}_mature.fastq) / 4 ))
    TOTAL_READS=\$(grep -oP "^\\d+" ${meta.id}_mapping.log | head -1)
    if [ -n "\$TOTAL_READS" ] && [ "\$TOTAL_READS" -gt 0 ]; then
        RATE=\$(python3 -c "print(f'{100 * \${TRNA_READS} / \${TOTAL_READS}:.1f}')")
        if python3 -c "exit(0 if \${TRNA_READS} / \${TOTAL_READS} < 0.02 else 1)"; then
            echo "" >&2
            echo "WARNING [${meta.id}]: Only \${TRNA_READS}/\${TOTAL_READS} reads (\${RATE}%) identified as tRNA." >&2
            echo "WARNING [${meta.id}]: This usually indicates incorrect adapter trimming." >&2
            echo "WARNING [${meta.id}]: Run 'detect_adapter.py <your.fastq.gz>' to check." >&2
            echo "" >&2
        fi
    fi

    gzip ${meta.id}_mature.fastq

    # Clean up intermediate SAM files
    rm -f ${meta.id}_artificial.sam ${meta.id}_filtered.sam
    """
}

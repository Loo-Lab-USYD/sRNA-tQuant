/*
 * CALL_MODIFICATIONS: variant calling on tRNA alignments to detect
 * post-transcriptional modifications, then adjust RPM if modifications
 * fall within the anticodon triplet.
 *
 * Only runs when params.skip_modification_calling is false.
 */

process CALL_VARIANTS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(sam)
    path cluster_fasta
    path cluster_fai

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: vcf

    script:
    """
    echo "call_variants v2" > /dev/null

    # Convert to sorted BAM with read groups
    samtools view -@ ${task.cpus} -bS ${sam} | samtools sort -@ ${task.cpus} -o ${meta.id}.bam
    samtools index -@ ${task.cpus} ${meta.id}.bam

    # Add read group (required by bcftools)
    samtools addreplacerg -@ ${task.cpus} -r '@RG\\tID:${meta.id}\\tSM:${meta.id}' \\
        -o ${meta.id}_rg.bam ${meta.id}.bam
    samtools index -@ ${task.cpus} ${meta.id}_rg.bam

    # Call variants
    bcftools mpileup --threads ${task.cpus} -f ${cluster_fasta} ${meta.id}_rg.bam | \\
        bcftools call --threads ${task.cpus} -mv --ploidy 1 -Oz -o ${meta.id}.vcf.gz

    bcftools index ${meta.id}.vcf.gz

    rm -f ${meta.id}.bam ${meta.id}.bam.bai ${meta.id}_rg.bam ${meta.id}_rg.bam.bai
    """
}

process COUNT_ALLELES {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(sam)
    tuple val(meta2), path(vcf)
    path cluster_fasta
    path cluster_fai

    output:
    tuple val(meta), path("${meta.id}.ASE.csv"), emit: ase

    script:
    """
    echo "count_alleles v2" > /dev/null

    # Convert to sorted BAM (using unique-only reads for mod calling)
    samtools view -@ ${task.cpus} -bS ${sam} | samtools sort -@ ${task.cpus} -o ${meta.id}.bam
    samtools index -@ ${task.cpus} ${meta.id}.bam

    # Count alleles at variant sites
    count_alleles.py ${meta.id}.bam ${vcf} ${cluster_fasta} ${meta.id}.ASE.csv --threads ${task.cpus}

    rm -f ${meta.id}.bam ${meta.id}.bam.bai
    """
}

process ADJUST_MODIFICATIONS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(rpm_csv)
    tuple val(meta2), path(ase_csv)
    path cluster_info
    path trna_list
    tuple val(meta3), path(cluster_rpm)

    output:
    tuple val(meta), path("${meta.id}_adjusted.RPM.csv"), emit: rpm

    script:
    def cluster_arg = cluster_rpm.size() > 0 ? "--cluster-rpm ${cluster_rpm}" : ''
    """
    echo "adjust_mods v2" > /dev/null
    adjust_modifications.py \\
        ${rpm_csv} \\
        ${ase_csv} \\
        ${cluster_info} \\
        ${trna_list} \\
        ${meta.id}_adjusted.RPM.csv \\
        ${cluster_arg}
    """
}

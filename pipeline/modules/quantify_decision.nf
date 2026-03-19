/*
 * QUANTIFY_DECISION: anticodon-aware multimapper resolution and RPM calculation.
 * Legacy-equivalent quantification path.
 */

process RESOLVE_MULTIMAPPERS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(sam)
    path cluster_info

    output:
    tuple val(meta), path("${meta.id}_resolved.sam"),      emit: sam
    tuple val(meta), path("${meta.id}_multimappers.sam"),   emit: multimappers

    script:
    """
    echo "resolve_multimappers v2" > /dev/null

    # Name-sort the SAM (required for grouping by read name)
    samtools sort -@ ${task.cpus} -n -O sam ${sam} > ${meta.id}_nsorted.sam

    # Resolve multimappers by anticodon identity
    resolve_multimappers.py \\
        ${meta.id}_nsorted.sam \\
        ${cluster_info} \\
        ${meta.id}_resolved.sam \\
        ${meta.id}_multimappers.sam

    rm -f ${meta.id}_nsorted.sam
    """
}

process QUANTIFY_RPM {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(sam)
    path cluster_info
    path trna_list

    output:
    tuple val(meta), path("${meta.id}.RPM.csv"),         emit: rpm
    tuple val(meta), path("${meta.id}.cluster_rpm.csv"), emit: cluster_rpm

    script:
    """
    echo "quantify_rpm v4" > /dev/null

    # Convert to BAM and get expression counts
    samtools view -@ ${task.cpus} -bS ${sam} | samtools sort -@ ${task.cpus} -o ${meta.id}.bam
    samtools index -@ ${task.cpus} ${meta.id}.bam
    samtools idxstats ${meta.id}.bam > ${meta.id}_expression.txt

    # Create empty ASE file (modification adjustment happens separately)
    echo -e "contig\\tposition\\tvariantID\\trefAllele\\taltAllele\\trefCount\\taltCount\\ttotalCount\\tlowMAPQDepth\\tlowBaseQDepth\\trawDepth\\totherBases\\timproperPairs" > ${meta.id}_empty.ASE.csv

    # Compute RPM per anticodon
    compute_expression.py \\
        ${meta.id}_expression.txt \\
        ${meta.id}_empty.ASE.csv \\
        ${cluster_info} \\
        ${trna_list} \\
        ${meta.id}.RPM.csv \\
        ${params.rpm_clamp_threshold ? "--clamp-threshold ${params.rpm_clamp_threshold}" : ''}

    # Empty cluster RPM file (Decision path doesn't have per-cluster Salmon RPMs)
    echo "cluster,rpm" > ${meta.id}.cluster_rpm.csv

    rm -f ${meta.id}.bam ${meta.id}.bam.bai
    """
}

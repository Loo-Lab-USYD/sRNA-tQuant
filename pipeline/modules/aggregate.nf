/*
 * AGGREGATE: merge per-sample RPM CSVs into a single matrix and run MultiQC.
 */

process AGGREGATE_RESULTS {
    tag "aggregate"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path rpm_csvs
    path ase_csvs

    output:
    path "results_rpm.csv",       emit: rpm_matrix
    path "results_rpm_nomod.csv", emit: rpm_nomod_matrix
    path "modifications.csv",     emit: modifications, optional: true

    script:
    def ase_arg = ase_csvs.name != 'NO_ASE' ? "--ase ${ase_csvs}" : ''
    """
    aggregate_results.py . ${rpm_csvs} ${ase_arg}
    """
}

process MULTIQC {
    tag "multiqc"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path multiqc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    """
    multiqc . --force --outdir .
    """
}

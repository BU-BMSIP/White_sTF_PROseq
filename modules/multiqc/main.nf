// File: modules/multiqc/multiqc.nf
process MULTIQC {
    label 'process_single'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    container 'ghcr.io/bf528/multiqc:latest'

    input:
    path fastqc_reports

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data

    script:
    """
    multiqc ${fastqc_reports} --outdir . --force
    """
}

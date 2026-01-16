process REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc --no-ai .
    """
}
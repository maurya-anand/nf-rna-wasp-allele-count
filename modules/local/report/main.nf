process REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc --no-ai .
    """
}
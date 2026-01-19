process ADAPTER_TRIM {
    publishDir "${params.outdir}/${meta.sampleid}/fastqc", mode: 'copy', pattern: "*_fastqc.{zip,html}"
    publishDir "${params.outdir}/${meta.sampleid}/trim_galore", mode: 'copy', pattern: "*trimming_report.txt"


    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    tuple val(meta), path("${meta.sampleid}_1_trimmed.fq.gz"), path("${meta.sampleid}_2_trimmed.fq.gz"), emit: reads
    path "*_fastqc.{zip,html}", emit: fastqc
    path("*_trimming_report.txt"), emit: log

    script:
    def fq1_base = fastq_1.toString().tokenize('.')[0]
    def fq2_base = fastq_2.toString().tokenize('.')[0]
    """
    total_threads=${task.cpus}
    set -euo pipefail
    trim_galore \
        --paired \
        --illumina \
        --fastqc \
        --cores \${total_threads} \
        ${fastq_1} ${fastq_2}
    mv ${fq1_base}_val_1.fq.gz ${meta.sampleid}_1_trimmed.fq.gz
    mv ${fq2_base}_val_2.fq.gz ${meta.sampleid}_2_trimmed.fq.gz
    mv ${fq1_base}_val_1_fastqc.html ${meta.sampleid}_1_val_1_fastqc.html
    mv ${fq1_base}_val_1_fastqc.zip ${meta.sampleid}_1_val_1_fastqc.zip
    mv ${fq2_base}_val_2_fastqc.html ${meta.sampleid}_2_val_2_fastqc.html
    mv ${fq2_base}_val_2_fastqc.zip ${meta.sampleid}_2_val_2_fastqc.zip
    mv ${fastq_1}_trimming_report.txt ${meta.sampleid}_1_trimming_report.txt
    mv ${fastq_2}_trimming_report.txt ${meta.sampleid}_2_trimming_report.txt
    """
}
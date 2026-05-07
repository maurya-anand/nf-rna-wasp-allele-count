#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { ADAPTER_TRIM } from './modules/local/adapter_trim'
include { STAR_GENOME_INDEX } from './modules/local/index_genome'
include { SUBSET_1KGP_VCF } from './modules/local/subset_1kgp'
include { STAR_ALIGNMENT_WASP } from './modules/local/alignment'
include { ALLELE_COUNT } from './modules/local/allele_count'
include { REPORT } from './modules/local/report'
include { CHECK_FASTQ } from './modules/local/check_fastq'

workflow {
    star_idx_ch = STAR_GENOME_INDEX(
        channel.fromPath(params.reference_fa),
        channel.fromPath(params.gencode_gtf),
    )
    reads_ch = channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def meta = [
                sampleid: row.sample
            ]
            [meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
        }
        .groupTuple(by: 0)

    vcf_ch_in = channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def meta = [
                sampleid: row.sample
            ]
            [meta, file(params.phased_vcf_dir)]
        }
        .unique()

    def regions_vcf_ch = params.regions_vcf ? channel.fromPath(params.regions_vcf, checkIfExists: true) : channel.empty()
    merged_reads_ch = CHECK_FASTQ(reads_ch)
    trimmed_reads_ch = ADAPTER_TRIM(merged_reads_ch)
    vcf_ch = SUBSET_1KGP_VCF(vcf_ch_in)
    align_in_ch = vcf_ch.subset_phased_vcf
        .join(trimmed_reads_ch.reads)
        .combine(star_idx_ch.star_index_dir)
        .map { meta, vcf, fq1, fq2, star_dir ->
            [meta, vcf, star_dir, fq1, fq2]
        }
    ac_in_ch = STAR_ALIGNMENT_WASP(align_in_ch)
    if (params.regions_vcf) {
        ALLELE_COUNT(
            ac_in_ch.bam,
            file(params.reference_fa),
            regions_vcf_ch,
        )
    }
    report_in_ch = channel.empty()
    report_in_ch = report_in_ch.mix(trimmed_reads_ch.fastqc)
    report_in_ch = report_in_ch.mix(trimmed_reads_ch.log.flatten())
    report_in_ch = report_in_ch.mix(ac_in_ch.log.map { _meta, log_file -> log_file })
    report_in_ch = report_in_ch.mix(ac_in_ch.stats.map { _meta, log_file -> log_file })
    REPORT(report_in_ch.collect())
}

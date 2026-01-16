#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOME_INDEX } from './modules/local/index_genome'
include { SUBSET_1KGP_VCF } from './modules/local/subset_1kgp'
include { STAR_ALIGNMENT_WASP } from './modules/local/alignment'
include { ALLELE_COUNT } from './modules/local/allele_count'
include { REPORT } from './modules/local/report'
workflow {
    star_idx_ch = STAR_GENOME_INDEX(
        channel.fromPath(params.reference_fa),
        channel.fromPath(params.gencode_gtf)
    )
    reads_ch = channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def meta = [
                sampleid: row.sample,
                fastq_1: file(row.fastq_1, checkIfExists: true),
                fastq_2: file(row.fastq_2, checkIfExists: true)
            ]
            [ meta, file(params.phased_vcf_dir) ]
        }
    vcf_ch = SUBSET_1KGP_VCF(reads_ch)
    reads_ch_rekeyed = reads_ch.map { meta, vcf_dir -> [ meta.sampleid, meta, vcf_dir ] }
    align_in_ch = vcf_ch.subset_phased_vcf
        .join(reads_ch_rekeyed)
        .combine(star_idx_ch.star_index_dir)
        .map { _sampleid, vcf, meta, _vcf_dir, star_dir ->
            [ meta, vcf, star_dir, meta.fastq_1, meta.fastq_2 ]
        }
    ac_in_ch = STAR_ALIGNMENT_WASP(align_in_ch)
    ALLELE_COUNT(
        ac_in_ch.bam,
        channel.fromPath(params.reference_fa),
        channel.fromPath(params.regions_vcf)
    )
    report_in_ch = channel.empty()
    report_in_ch = report_in_ch.mix(ac_in_ch.log.map { _meta, log_file -> log_file })
    report_in_ch = report_in_ch.mix(ac_in_ch.stats.map { _meta, log_file -> log_file })
    REPORT(report_in_ch.collect())
}
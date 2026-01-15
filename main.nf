#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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
            [ meta ]
        }
    vcf_input_ch = reads_ch
        .map { meta -> [ meta, file(params.phased_vcf_dir) ] }
    vcf_ch = SUBSET_1KGP_VCF(vcf_input_ch)
    align_in_ch = vcf_ch.subset_phased_vcf
        .join(reads_ch)
        .combine(star_idx_ch.star_index_dir)
    STAR_ALIGNMENT_WASP(align_in_ch)
}

process STAR_ALIGNMENT_WASP {
    input:
    tuple val(meta), path(phased_1KGP_vcf), path(phased_1KGP_vcf_idx), path(star_index_dir)
    output:
    tuple val(meta), path("${meta.sampleid}.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.sampleid}.Log.final.out"), emit: log
    tuple val(meta), path("${meta.sampleid}.SJ.out.tab"), emit: sj_out
    script:
    """
    set -euo pipefail

    STAR \\
    --runMode alignReads \\
    --runThreadN ${task.cpus} \\
    --twopassMode Basic \\
    --genomeDir ${star_index_dir} \\
    --varVCFfile ${phased_1KGP_vcf} \\
    --waspOutputMode SAMtag \\
    --readFilesIn ${meta.fastq_1} ${meta.fastq_2} \\
    --readFilesCommand zcat \\
    --outSAMtype BAM SortedByCoordinate \\
    --outSAMunmapped Within \\
    --outFileNamePrefix ${meta.sampleid}. \\
    --outFilterMultimapNmax 20 \\
    --alignSJoverhangMin 8 \\
    --alignSJDBoverhangMin 1 \\
    --outFilterMismatchNmax 999 \\
    --outFilterMismatchNoverLmax 0.1 \\
    --alignIntronMin 20 \\
    --alignIntronMax 1000000 \\
    --alignMatesGapMax 1000000 \\
    --outFilterType BySJout \\
    --outFilterScoreMinOverLread 0.33 \\
    --outFilterMatchNminOverLread 0.33 \\
    --outFilterMatchNmin 0 \\
    --limitSjdbInsertNsj 1200000 \\
    --outSAMstrandField intronMotif \\
    --outFilterIntronMotifs None \\
    --alignSoftClipAtReferenceEnds Yes \\
    --quantMode TranscriptomeSAM GeneCounts \\
    --outSAMattrRGline ID:${meta.sampleid} SM:${meta.sampleid} \\
    --outSAMattributes NH HI AS nM NM ch \\
    --chimOutJunctionFormat 0 \\
    --chimSegmentMin 15 \\
    --chimJunctionOverhangMin 15 \\
    --chimOutType Junctions WithinBAM SoftClip \\
    --chimMainSegmentMultNmax 1 \\
    --genomeLoad NoSharedMemory
    """
}

process STAR_GENOME_INDEX {
    input:
    path ref_fa
    path gtf

    output:
    path("STAR_INDEX"), emit: star_index_dir

    script:
    """
    set -euo pipefail

    mkdir -p STAR_INDEX
    STAR --runMode genomeGenerate \\
      --runThreadN ${task.cpus} \\
      --genomeDir STAR_INDEX \\
      --genomeFastaFiles ${ref_fa} \\
      --sjdbGTFfile ${gtf}
    """
}

process SUBSET_1KGP_VCF {
    input:
    tuple val(meta), path(phased_vcf_dir)
    output:
    tuple val(meta.sampleid), path("${meta.sampleid}.1KGP.snps.het.vcf.gz"), path("${meta.sampleid}.1KGP.snps.het.vcf.gz.tbi"), emit: subset_phased_vcf
    script:
    """
    total_threads=${task.cpus}
    if ls ${phased_vcf_dir}/*.vcf.gz 1> /dev/null 2>&1; then
        echo "Subsetting from Phased VCF dir"
        for i in {1..22}; do
            chrom=chr\${i}
            vcf_file=${phased_vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_\${chrom}.filtered.shapeit2-duohmm-phased.vcf.gz
            bcftools view --threads \${total_threads} --no-update -s ${meta.sampleid} -v snps \${vcf_file} | bcftools view --threads \${total_threads} --no-update -e 'GT=".|."' -Oz -o ${meta.sampleid}.\${chrom}.snps.vcf.gz
            bcftools view --threads \${total_threads} --no-update -i 'GT="het"' ${meta.sampleid}.\${chrom}.snps.vcf.gz | bcftools norm -m+ | bcftools view --threads \${total_threads} -m2 -M2 -Oz -o ${meta.sampleid}.\${chrom}.snps.het.vcf.gz
            tabix -p vcf ${meta.sampleid}.\${chrom}.snps.het.vcf.gz
        done

        chrom=chrX
        vcf_file=${phased_vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_\${chrom}.filtered.eagle2-phased.v2.vcf.gz
        bcftools view --threads \${total_threads} --no-update -s ${meta.sampleid} -v snps \${vcf_file} | bcftools view --threads \${total_threads} --no-update -e 'GT=".|."' -Oz -o ${meta.sampleid}.\${chrom}.snps.vcf.gz
        bcftools view --threads \${total_threads} --no-update -i 'GT="het"' ${meta.sampleid}.\${chrom}.snps.vcf.gz | bcftools norm -m+ | bcftools view --threads \${total_threads} -m2 -M2 -Oz -o ${meta.sampleid}.\${chrom}.snps.het.vcf.gz
        tabix -p vcf ${meta.sampleid}.\${chrom}.snps.het.vcf.gz
        {
            for i in \$(seq 1 22); do echo "${meta.sampleid}.chr\${i}.snps.het.vcf.gz"; done
            [[ -s "${meta.sampleid}.chrX.snps.het.vcf.gz" ]] && echo "${meta.sampleid}.chrX.snps.het.vcf.gz"
        } > list.txt
        bcftools concat -f list.txt -Ou \\
        | bcftools sort -Oz -o ${meta.sampleid}.1KGP.snps.het.vcf.gz

        tabix -p vcf ${meta.sampleid}.1KGP.snps.het.vcf.gz
    else
        echo "Phased VCF dir is empty. Creating empty output files."
        touch ${meta.sampleid}.1KGP.snps.het.vcf.gz
        touch ${meta.sampleid}.1KGP.snps.het.vcf.gz.tbi
    fi
    """
}

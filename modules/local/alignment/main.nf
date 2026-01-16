process STAR_ALIGNMENT_WASP {
    publishDir "${params.outdir}/${meta.sampleid}/alignment", mode: 'copy'

    input:
    tuple val(meta), path(phased_1KGP_vcf), path(star_index_dir), path(fastq_1), path(fastq_2)

    output:
    tuple val(meta), path("${meta.sampleid}.Aligned.sortedByCoord.out.bam"), path("${meta.sampleid}.Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(meta), path("${meta.sampleid}.samtools_stats.txt"), emit: stats
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
    --readFilesIn ${fastq_1} ${fastq_2} \\
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

    samtools index ${meta.sampleid}.Aligned.sortedByCoord.out.bam
    samtools stats ${meta.sampleid}.Aligned.sortedByCoord.out.bam > ${meta.sampleid}.samtools_stats.txt
    """
}
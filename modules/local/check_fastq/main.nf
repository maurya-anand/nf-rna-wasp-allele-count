process CHECK_FASTQ {
    tag "$meta.sampleid"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("*_1.fastq.gz"), path("*_2.fastq.gz")

    script:
    if (reads1 instanceof List && reads1.size() > 1) {
        """
        cat ${reads1.join(' ')} > ${meta.sampleid}_1.fastq.gz
        cat ${reads2.join(' ')} > ${meta.sampleid}_2.fastq.gz
        """
    } else {
        def r1 = reads1 instanceof List ? reads1[0] : reads1
        def r2 = reads2 instanceof List ? reads2[0] : reads2
        """
        ln -s ${r1} ${meta.sampleid}_1.fastq.gz
        ln -s ${r2} ${meta.sampleid}_2.fastq.gz
        """
    }
}
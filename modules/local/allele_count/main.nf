process ALLELE_COUNT{
    publishDir "${params.outdir}/${meta.sampleid}/allele_counts", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference_fa
    path regions_vcf

    output:
    tuple val(meta.sampleid),
        path("${meta.sampleid}.allele_counts.with_qual.tsv"),
        path("${meta.sampleid}.allele_counts.with_qual.tsv.gz"),
        path("${meta.sampleid}.allele_counts.with_qual.tsv.gz.tbi"),
        emit: allele_counts

    script:
    """
    set -euo pipefail

    [[ -f "${reference_fa}.fai" ]] || samtools faidx ${reference_fa}

    bgzip -c ${regions_vcf} > sites.vcf.gz

    tabix -p vcf sites.vcf.gz

    echo -e "CHROM\\tPOS\\tREF\\tALT\\tDP\\tREF_COUNT\\tALT_COUNT\\tVAF\\tMQ\\tQS\\tQSsum" \
        > ${meta.sampleid}.allele_counts.with_qual_and_vaf.tsv
    
    bcftools mpileup \\
        -f ${reference_fa} \\
        -R sites.vcf.gz \\
        -a FORMAT/AD,FORMAT/DP,INFO/MQ,INFO/QS,INFO/QSsum \\
        -Q 20 -q 20 \\
        -Ou ${bam} \\
    | bcftools call -m none -Ou \\
    | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%DP]\\t[%AD]\\t%MQ\\t%QS\\t%QSsum\\n' \\
    | awk -F '\\t' 'BEGIN{OFS="\\t"}
        {
        split(\$6, ad, ",");
        ref=ad[1]; alt=ad[2];
        vaf = (ref+alt > 0) ? alt/(ref+alt) : "NA";
        print \$1,\$2,\$3,\$4,\$5,ref,alt,vaf,\$7,\$8,\$9
        }' \
    >> ${meta.sampleid}.allele_counts.with_qual_and_vaf.tsv

    bgzip -c ${meta.sampleid}.allele_counts.with_qual_and_vaf.tsv \
        > ${meta.sampleid}.allele_counts.with_qual_and_vaf.tsv.gz

    tabix -s 1 -b 2 -e 2 \
        ${meta.sampleid}.allele_counts.with_qual_and_vaf.tsv.gz
    """
}
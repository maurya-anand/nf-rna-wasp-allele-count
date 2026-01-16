process ALLELE_COUNT{
    publishDir "${params.outdir}/${meta.sampleid}/allele_counts", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference_fa
    path regions_vcf

    output:
    tuple val(meta.sampleid), path("${meta.sampleid}.vcf"), path("${meta.sampleid}.allele_counts_vaf.tsv"), emit: allele_counts

    script:
    """
    set -euo pipefail

    [[ -f "${reference_fa}.fai" ]] || samtools faidx ${reference_fa}

    bgzip -c ${regions_vcf} > sites.vcf.gz

    tabix -p vcf sites.vcf.gz

    echo -e "CHROM\\tPOS\\tREF\\tALT\\tDP\\tREF_COUNT\\tALT_COUNT\\tVAF" \
        > ${meta.sampleid}.allele_counts_vaf.tsv
    
    bcftools mpileup \\
        -f ${reference_fa} \\
        -R sites.vcf.gz \\
        -a FORMAT/AD,FORMAT/DP \\
        -Ou ${bam} \\
    | bcftools call -m -Ou -o ${meta.sampleid}.vcf
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%DP]\\t[%AD]\\n' ${meta.sampleid}.vcf \\
    | awk -F '\\t' 'BEGIN{OFS="\\t"}
        {
        split(\$6, ad, ",");
        ref=ad[1]; alt=ad[2];
        vaf = (ref+alt > 0) ? alt/(ref+alt) : "NA";
        print \$1,\$2,\$3,\$4,\$5,ref,alt,vaf
        }' \
    >> ${meta.sampleid}.allele_counts_vaf.tsv
    """
}
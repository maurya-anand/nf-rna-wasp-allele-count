process SUBSET_1KGP_VCF {
    publishDir "${params.outdir}/${meta.sampleid}/extracted_1KGP_het_snps", mode: 'copy'

    input:
    tuple val(meta), path(phased_vcf_dir)

    output:
    tuple val(meta), path("${meta.sampleid}.1KGP.snps.het.vcf"), emit: subset_phased_vcf

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
        | bcftools sort --temp-dir ./ -Ov -o ${meta.sampleid}.1KGP.snps.het.vcf

    else
        echo "Phased VCF dir is empty. Creating empty output files."
        touch ${meta.sampleid}.1KGP.snps.het.vcf
    fi
    """
}
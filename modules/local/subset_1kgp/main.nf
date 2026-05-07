process SUBSET_1KGP_VCF {
    publishDir "${params.outdir}/${meta.sampleid}/extracted_1KGP_het_snps", mode: 'copy'

    input:
    tuple val(meta), path(phased_vcf_dir)

    output:
    tuple val(meta), path("${meta.sampleid}.1KGP.snps.het.vcf"), emit: subset_phased_vcf

    script:
    """
    export TMPDIR=${task.workDir}

    total_threads=${task.cpus}

    threads_per_job=\$(( total_threads > 2 ? 2 : total_threads ))
    parallel_jobs=\$(( total_threads / threads_per_job ))
    parallel_jobs=\$(( parallel_jobs < 1 ? 1 : parallel_jobs ))

    process_chrom() {
        local chrom=\$1
        local vcf_file=\$2
        local sample=\$3
        local threads=\$4

        bcftools view --threads \${threads} --no-update -s \${sample} -v snps \${vcf_file} \\
            | bcftools view --threads \${threads} --no-update -e 'GT=".|."' -Oz -o \${sample}.\${chrom}.snps.vcf.gz

        bcftools view --threads \${threads} --no-update -i 'GT="het"' \${sample}.\${chrom}.snps.vcf.gz \\
            | bcftools norm -m+ \\
            | bcftools view --threads \${threads} -m2 -M2 -Oz -o \${sample}.\${chrom}.snps.het.vcf.gz

        tabix -p vcf \${sample}.\${chrom}.snps.het.vcf.gz
    }
    export -f process_chrom

    if ls ${phased_vcf_dir}/*.vcf.gz 1> /dev/null 2>&1; then
        echo "Subsetting from Phased VCF dir"

        seq 1 22 | xargs -P \${parallel_jobs} -I{} bash -c '
            process_chrom \\
                "chr{}" \\
                "${phased_vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{}.filtered.shapeit2-duohmm-phased.vcf.gz" \\
                "${meta.sampleid}" \\
                "\${threads_per_job}"
        '
        process_chrom \\
            "chrX" \\
            "${phased_vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz" \\
            "${meta.sampleid}" \\
            "\${total_threads}"

        {
            for i in \$(seq 1 22); do echo "${meta.sampleid}.chr\${i}.snps.het.vcf.gz"; done
            [[ -s "${meta.sampleid}.chrX.snps.het.vcf.gz" ]] && echo "${meta.sampleid}.chrX.snps.het.vcf.gz"
        } > list.txt

        bcftools concat --threads \${total_threads} -f list.txt -Ou \\
            | bcftools sort --temp-dir ${task.workDir} -Ov -o ${meta.sampleid}.1KGP.snps.het.vcf
    else
        echo "Phased VCF dir is empty. Creating empty output files."
        touch ${meta.sampleid}.1KGP.snps.het.vcf
    fi
    """
}

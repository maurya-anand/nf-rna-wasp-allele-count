# nf-rna-wasp-allele-count

A Nextflow DSL2 pipeline to perform allele-specific expression analysis from RNA-seq data using a WASP-aware alignment strategy.

> [!NOTE]
> This pipeline adopts the genotype-aware alignment strategy used in:
>
> Taylor DJ, Chhetri SB, Tassia MG, et al. *Sources of gene expression variation in a globally diverse human cohort.* **Nature** 632, 122–130 (2024). <https://doi.org/10.1038/s41586-024-07708-2>
>
> This workflow uses STAR alignment with WASP correction to prevent reference-mapping bias at heterozygous sites. While the original study focused on regulatory and splicing variation, this pipeline extends the same alignment principles to explicit allele-specific read counting from RNA-seq data.

## Features

- Genotype-aware RNA-seq alignment with **STAR + WASP correction**.
- Allele-specific read counting (REF/ALT) from RNA-seq BAMs.
- Reports total depth and variant allele fraction.

## Requirements

- Nextflow **≥ 22.10.0**
- Singularity / Apptainer or Docker (or compatible container runtime)
- Phased VCFs for each sample (e.g., 1000G high-coverage release)
- Reference genome FASTA + GTF matching alignment build

## Quick Start

### 1. Prepare your sample sheet (CSV)

`samplesheet.csv`

```csv
sample,fastq_1,fastq_2
sampleID,sampleID_1.fastq.gz,sampleID_2.fastq.gz
```

### 2. Download reference data

Download the high-coverage phased WGS VCFs an store all the files in a directory:

<https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased>

Expected output:

```bash
1000G_2504_high_coverage_20201028_3202_phased_vcfs/
├── CCDG_14151_B01_GRM_WGS_2020-08-05_chr*.filtered.shapeit2-duohmm-phased.vcf.gz
├── CCDG_14151_B01_GRM_WGS_2020-08-05_chr*.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
├── CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz
└── CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz.tbi
```

### 3. Reference genome and annotation

```bash
mkdir reference
cd reference
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gunzip gencode.v38.annotation.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

Expected output:

```bash
reference/
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
└── gencode.v38.annotation.gtf
```

### 4. Define SNPs of interest (VCF)

The pipeline expects a VCF so REF and ALT alleles are explicit.

`query_sites.vcf`

```tsv
##fileformat=VCFv4.2
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
chr1    123 .   C   T   .   .  .
```

### 5. Set pipeline parameters

```groovy
params {
    sample_sheet = samplesheet.csv
    phased_vcf_dir = 1000G_2504_high_coverage_20201028_3202_phased_vcfs
    reference_fa = GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    gencode_gtf = gencode.v38.annotation.gtf
    regions_vcf = query_sites.vcf
    outdir = results/
}
```

### 6. Run the pipeline

```bash
nextflow run main.nf -profile docker
```

## Workflow Overview

The pipeline consists of the following main steps:

- Genome index generation (STAR_GENOME_INDEX)
  - Builds a STAR genome index from GRCh38 and GENCODE v38.
  - Executed once and reused across all samples.
  - Output: `STAR_INDEX/`

- Extraction of phased variants (SUBSET_1KGP_VCF)
  - Subsets 1000 Genomes phased WGS VCFs to a single sample.
  - Extracts biallelic heterozygous SNPs only.
  - Output: `sampleID.1KGP.snps.het.vcf`

- Alignment with WASP correction (STAR_ALIGNMENT_WASP)
  - Aligns reads using STAR with genotype-aware variant input
  - Applies WASP correction to mitigate reference mapping bias
  - Uses per-sample phased variants via `--varVCFfile`.
  - Output: `sampleID.Aligned.sortedByCoord.out.bam` (coordinate sorted BAM file with WASP tags)

- Allele-specific read counting (ALLELE_COUNT)
  - Uses bcftools mpileup restricted to SNPs of interest.
  - Reports:
    - Total depth (DP)
    - REF and ALT read counts (AD)
    - Computes variant allele fraction (VAF = ALT / (REF + ALT)).
  - Output: `sampleID.allele_counts.with_qual_and_vaf.tsv.gz`

## Customization

Resource requirements, containers, and execution profiles can be adjusted in:

- `nextflow.config`
- `conf/base.config`

## Components

Tools:

| Component | Version |
|-----------|---------|
| STAR      | 2.7.11b |
| SAMTOOLS  | 1.21    |
| BCFTOOLS  | 1.23    |
| TABIX     | 1.11    |

Container Image:

- `community.wave.seqera.io/library/bcftools_samtools_star_tabix:e294cd9e3fb171ce`

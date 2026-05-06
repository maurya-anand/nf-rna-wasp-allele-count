# Changelog

All notable changes to this project will be documented in this file.

## [1.0.2] - 2026-05-06

- Fixed HPC profile
  - Added SLURM executor throttling (queueSize, submitRateLimit) for large-scale runs.
  - Added retry logic with resource scaling for OOM and timeout failures.
  - Added scratch support for node-local I/O performance.
  - Capped memory and CPU resources with Math.min() to prevent exceeding node limits.
  - Removed redundant executor and queue declarations from process selectors.

## [1.0.1] - 2026-03-30

- Conditional execution of alelle count process.
- Added pipeline execution report.

## [1.0.0] - 2026-01-19

### Release Notes

- Initial release of nf-rna-wasp-allele-count, a Nextflow pipeline to perform allele-specific expression analysis from RNA-seq data using a WASP-aware alignment strategy.
- Implements modular Nextflow DSL2 workflow with sample-wise output organization.
- Major features:
  - Genotype-aware RNA-seq alignment with **STAR + WASP correction**.
  - Allele-specific read counting (REF/ALT) from RNA-seq BAMs.
  - Reporting of depth, mapping quality, base-quality summaries, and VAF.
- Supports execution with Docker, Singularity, Apptainer, and SLURM HPC environments.

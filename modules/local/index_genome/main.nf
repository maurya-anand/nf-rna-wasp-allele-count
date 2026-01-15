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
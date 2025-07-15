process PROSEQ_METRICS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bed), val(tss), val(tes)

    output:
    tuple val(sample_id),
          path("${sample_id}.binned_counts.tsv"),
          path("${sample_id}.metrics.txt")

    conda """
      channels:
        - conda-forge
        - bioconda
      dependencies:
        - python=3.10
        - pandas
        - numpy
    """

    script:
    """
    python3 ${projectDir}/scripts/proseq_count.py \
            $bed $tss $tes ${params.bin_size} $sample_id
    """
}

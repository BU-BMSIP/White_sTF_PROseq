process UMI_EXTRACT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1)

    output:
    tuple val(sample_id), path("${sample_id}.umi.fastq.gz")

    conda "bioconda::umi_tools=1.1.5"

    script:
    """
    umi_tools extract \
          --stdin $r1 \
          --stdout ${sample_id}.umi.fastq.gz \
          --bc-pattern=${params.umi_pattern} \
          --filter-cell-barcode \
          --extract-method=string
    """
}

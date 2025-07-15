process BAM2BED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bed")

    conda "bioconda::bedtools=2.31.1"

    script:
    """
    bedtools bamtobed -i $bam > ${sample_id}.bed
    """
}

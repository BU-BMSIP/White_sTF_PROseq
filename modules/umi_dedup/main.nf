process UMI_DEDUP {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam")

    conda "bioconda::umi_tools=1.1.5"

    script:
    """
    umi_tools dedup -I $bam -S ${sample_id}.dedup.bam \
                    --paired False --log=${sample_id}.dedup.log
    """
}

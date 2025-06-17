process BAMCOVERAGE {
    label 'process_medium'
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    container 'ghcr.io/bf528/deeptools:latest'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}.bw"), emit: bigwig

    script:
    """
    bamCoverage \
        --bam ${bam} \
        --outFileName ${sample_id}.bw \
        --normalizeUsing CPM \
        --binSize 10 \
        --numberOfProcessors ${task.cpus}
    """
}

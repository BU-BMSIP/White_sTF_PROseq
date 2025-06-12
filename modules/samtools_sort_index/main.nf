process SAMTOOLS_SORT_INDEX {
    label 'process_low'
    publishDir "${params.outdir}/aligned", mode: 'copy'
    

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: sorted_bam

    script:
    """
    samtools sort -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """
}

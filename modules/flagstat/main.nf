process FLAGSTAT {
    label 'process_single'
    publishDir "${params.outdir}/qc/flagstat", mode: 'copy'

    input:
        tuple val(sample_id), path(bam)

    output:
        path("${sample_id}.flagstat.txt"), emit: stats

    script:
    """
    samtools flagstat ${bam} > ${sample_id}.flagstat.txt
    """
}

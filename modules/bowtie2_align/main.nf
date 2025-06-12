process BOWTIE2_ALIGN {
    label      'process_medium'
    publishDir "${params.outdir}/aligned", mode: 'copy'
    

    input:
        path index_files
        tuple val(sample_id), path(reads), val(index_base)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), emit: bam

    script:
    def nthreads = task.cpus - 1
    """
    bowtie2 -p ${nthreads} -x ${index_base} -U ${reads} | \
        samtools view -Sb - > ${sample_id}.bam
    """

}
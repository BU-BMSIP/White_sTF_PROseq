process BOWTIE2_INDEX {
    label       'process_low'
    publishDir  "${params.outdir}/bowtie2/index", mode:'copy'
    container   'ghcr.io/bf528/bowtie2:latest'

    input:
        path fasta

    output:
        path "${fasta.baseName}.*.bt2", emit: index_files
        val(fasta.baseName),           emit: index_base

    script:
    """
    bowtie2-build --threads ${task.cpus} ${fasta} ${fasta.baseName}
    """
}
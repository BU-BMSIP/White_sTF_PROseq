process CUTADAPT_TRIM {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz")

    conda "bioconda::cutadapt=4.4"

    script:
    """
    cutadapt   -a ${params.adapter_r1} \
               -m ${params.min_len} \
               -q ${params.qual_cutoff} \
               -o ${sample_id}.trimmed.fastq.gz \
               $reads
    """
}

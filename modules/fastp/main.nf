process FASTP {
    label 'process_low'
    publishDir "${params.outdir}/fastp", mode: 'copy'
    container 'quay.io/biocontainers/fastp:0.23.2--h5f740d0_3'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.trimmed.fastq.gz"), emit: trimmed
    path("fastp.json"), emit: json
    path("fastp.html"), emit: html

    script:
    def out_prefix = sample_id + ".trimmed"
    """
    fastp \
        -i ${reads} \
        -o ${out_prefix}.fastq.gz \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --detect_adapter_for_pe \
        --thread 4 \
        --html fastp.html \
        --json fastp.json \
        --qualified_quality_phred 15 \
        --length_required 20
    """
}

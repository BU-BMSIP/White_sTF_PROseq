#!/usr/bin/env nextflow

process FASTQC {
   label 'process_single'
   publishDir "${params.outdir}/fastqc", mode: 'copy'
   container 'ghcr.io/bf528/fastqc:latest'
   publishDir params.outdir, mode: 'copy'
  
   input:
   tuple val(sample_id), path(fastq)

   output:
   path("*_fastqc.zip"), emit: zip
   path("*_fastqc.html"), emit: html
   

   script:
   """
   fastqc ${fastq}
   """
}
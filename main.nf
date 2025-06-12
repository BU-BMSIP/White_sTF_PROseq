#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ------------------------------------------------------------
// Parameter definitions
// ------------------------------------------------------------
params.samplesheet = "${params.samplesheet ?: 'samples.csv'}"
params.outdir      = "${params.outdir ?: 'results'}"
params.skip_index  = params.skip_index ?: false   // flag to reuse existing .bt2 files

// ------------------------------------------------------------
// Module imports
// ------------------------------------------------------------
include { FASTQC }        from './modules/fastqc/main.nf'
include { FASTP }         from './modules/fastp/main.nf'
include { MULTIQC }       from './modules/multiqc/main.nf'
include { BOWTIE2_INDEX } from './modules/bowtie2_index/main.nf'
include { BOWTIE2_ALIGN } from './modules/bowtie2_align/main.nf'
include { SAMTOOLS_SORT_INDEX } from './modules/samtools_sort_index/main.nf'
include { FLAGSTAT } from './modules/flagstat/main.nf'



// ------------------------------------------------------------
// Workflow definition
// ------------------------------------------------------------
workflow {

    //------------------------------------------------------------------
    // 1) Ingest sample sheet -> Channel of (sample_id, fastq_path)
    //------------------------------------------------------------------
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map   { row -> tuple( row.name, file(row.path) ) }
        .set   { fastq_ch }

    //------------------------------------------------------------------
    // 2) Raw-read QC
    //------------------------------------------------------------------
    def fastqc_out  = FASTQC( fastq_ch )
    def fastqc_zips = fastqc_out.zip.collect()

    MULTIQC( fastqc_zips )

    //------------------------------------------------------------------
    // 3) Read trimming
    //------------------------------------------------------------------
    def fastp_out     = FASTP( fastq_ch )
    def trimmed_reads = fastp_out.trimmed  // (sample_id, trimmed_fastq)

    //------------------------------------------------------------------
    // 4) Build (or reuse) Bowtie2 index
    //------------------------------------------------------------------
    def index_base
    def index_files

    if( params.skip_index ) {
        log.info "↪︎  Skipping index build – using pre‑existing merged_reference.*.bt2 files"
        index_base  = Channel.value('refs/merged_reference')
        index_files = Channel.fromPath('refs/merged_reference*.bt2')
    }
    else {
        def bowtie_idx  = BOWTIE2_INDEX( file('refs/merged_reference.fa') )
        index_files = bowtie_idx.index_files
        index_base  = bowtie_idx.index_base
    }

    //------------------------------------------------------------------
    // 5) Alignment
    //------------------------------------------------------------------
    def align_input   = trimmed_reads.combine(index_base)  // (sample_id, fastq, index_base)
    def aligned_bams  = BOWTIE2_ALIGN( index_files, align_input )

    // 6) Wrap your aligned BAMs:
    aligned_bams.map { sample_id, bam -> tuple(sample_id, bam) } \
        .set { unsorted_bams }

    sorted_bams = SAMTOOLS_SORT_INDEX(unsorted_bams)

    // 7) Get stats on sorted BAMs:
    // Only pass (sample_id, sorted.bam) to FLAGSTAT
    sorted_bams.map { sample_id, bam, bai -> tuple(sample_id, bam) } \
        .set { sorted_bams_for_flagstat }

    flagstat_output = FLAGSTAT(sorted_bams_for_flagstat)


}
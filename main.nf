#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


/* ─────────────────────────  PARAMETERS  ───────────────────────── */
params.samplesheet   = params.samplesheet ?: 'samples.csv'
params.outdir        = params.outdir      ?: 'results'
params.skip_index    = params.skip_index  ?: false

// dm6 spike-in
params.dm6_index     = params.dm6_index   ?: 'refs/dm6/dm6'   // basename *without* .bt2

// PRO-seq-specific
params.adapter_r1    = params.adapter_r1  ?: 'TGGAATTCTCGG'
params.min_len       = params.min_len     ?: 26
params.qual_cutoff   = params.qual_cutoff ?: 10
params.umi_pattern   = params.umi_pattern ?: 'NNNNNN'
params.skip_dedup    = params.skip_dedup  ?: false
params.bin_size      = params.bin_size    ?: 10
params.tss_bed       = params.tss_bed     ?: 'refs/plasmid_tss.bed'


/* ─────────────────────────  MODULE IMPORTS  ───────────────────── */
include { FASTQC }                 from './modules/fastqc/main.nf'
include { MULTIQC }                from './modules/multiqc/main.nf'

include { CUTADAPT_TRIM }          from './modules/cutadapt_trim/main.nf'
include { UMI_EXTRACT }            from './modules/umi_extract/main.nf'
include { UMI_DEDUP }              from './modules/umi_dedup/main.nf'

include { BOWTIE2_INDEX }          from './modules/bowtie2_index/main.nf'
include { BOWTIE2_ALIGN }          from './modules/bowtie2_align/main.nf'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_DM6 } from './modules/bowtie2_align/main.nf'

include { SAMTOOLS_SORT_INDEX }    from './modules/samtools_sort_index/main.nf'
include { FLAGSTAT }               from './modules/flagstat/main.nf'
include { FLAGSTAT as FLAGSTAT_DM6 }from './modules/flagstat/main.nf'

include { BAM2BED }                from './modules/bam2bed/main.nf'
include { PROSEQ_METRICS }         from './modules/proseq_metrics/main.nf'

include { BAMCOVERAGE }            from './modules/bamcoverage/main.nf'
include { BIGWIG_CORRELATION }     from './modules/bigwig_correlation/main.nf'



/* ───────────────────────────  WORKFLOW  ───────────────────────── */
workflow {

    /* 1 ─ Sample sheet  (sample_id , fastq path) */
    Channel
        .fromPath( params.samplesheet )
        .splitCsv(header:true)
        .map   { row -> tuple(row.name, file(row.path)) }
        .set   { fastq_ch }


    /* 2 ─ FastQC + MultiQC (run #1) */
    fastqc_out  = FASTQC( fastq_ch )

    fastqc_zip_list = fastqc_out.zip          // property is already a Channel
                         .collect()           // → single List[path]
    MULTIQC( fastqc_zip_list )


    /* 3 ─ Adapter trimming */
    trimmed_reads = CUTADAPT_TRIM( fastq_ch )         // (sample, trimmed.fq.gz)


    /* 4 ─ UMI extraction */
    umi_reads = UMI_EXTRACT( trimmed_reads )          // (sample, umi.fq.gz)


    /* 5 ─ Bowtie2 index (build or reuse) */
    def index_files
    def index_base
    if( params.skip_index ) {
        index_files = Channel.fromPath('refs/merged_reference*.bt2')
        index_base  = Channel.value      ('refs/merged_reference')
    }
    else {
        built       = BOWTIE2_INDEX( file('refs/merged_reference.fa') )
        index_files = built.index_files
        index_base  = built.index_base                      // value channel
    }


    /* 6 ─ Alignment to merged reference
            ‣ module expects:
              - channel 1: *.bt2 shards
              - channel 2: (sample, fq , basename)
    */
    align_input = umi_reads.combine(index_base)          // add basename
    aligned_bams = BOWTIE2_ALIGN( index_files, align_input )


    /* 6b ─ Alignment to dm6 spike-in */
    dm6_idx_files = Channel.fromPath("${params.dm6_index}*.bt2")
    dm6_idx_base  = Channel.value( params.dm6_index )
    dm6_input     = umi_reads.combine( dm6_idx_base )

    dm6_bams  = BOWTIE2_ALIGN_DM6( dm6_idx_files, dm6_input )
    FLAGSTAT_DM6( dm6_bams.map { id,bam -> tuple(id,bam) } )


    /* 7 ─ Sort + index primary BAM */
    sorted_bams = SAMTOOLS_SORT_INDEX(
                     aligned_bams.map { id,bam -> tuple(id,bam) }
                  )                            // (sample, sorted.bam, .bai)


    /* 8 ─ UMI de-dup (optional) */
    final_bams = params.skip_dedup \
        ? sorted_bams \
        : UMI_DEDUP( sorted_bams.map{ id,bam,bai -> tuple(id,bam) } )


    /* 9 ─ Flagstat on final BAM */
    FLAGSTAT( final_bams.map{ id,bam -> tuple(id,bam) } )


    /* 10 ─ bigWig coverage */
    bigwig_ch = BAMCOVERAGE( final_bams )       // (sample, *.bw)


    /* 11 ─ Single-nt trace + pausing index */
    bed_ch = BAM2BED( final_bams.map{ id,bam -> tuple(id,bam) } )

    coords_ch = Channel
                  .fromPath( params.tss_bed )
                  .splitCsv(sep:'\t', header:true)
                  .map { row -> tuple(row.name, row.tss as int, row.tes as int) }

    metrics_ch = PROSEQ_METRICS( bed_ch.join(coords_ch) )


    /* 12 ─ bigWig correlation */
    BIGWIG_CORRELATION( bigwig_ch.collect(), file('refs/genes.bed') )


    /* 13 ─ Final MultiQC (FastQC + pausing metrics) */
def mqc_inputs = fastqc_zip_list              \
                   .concat(                   \
                     metrics_ch.map{ s,tsv,txt -> txt } \
                   )
MULTIQC( mqc_inputs.collect() )

}

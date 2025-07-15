process BIGWIG_CORRELATION {

    label 'process_low'
    publishDir "${params.outdir}/correlation", mode: 'copy'

    input:
    path bigwigs, stageAs: 'bigwig_*'
    path bed_file

    output:
    path "correlation_matrix.tsv", emit: matrix
    path "correlation_heatmap.png", emit: heatmap
    path "summary_counts.npz",     emit: summary_npz

    script:
    """
    # deepTools v3+
    multiBigwigSummary BED-file \
        --bwfiles ${bigwigs.join(' ')} \
        --BED     ${bed_file} \
        --outFileName   summary_counts.npz \
        --outRawCounts  summary_counts.tsv \
        -p ${task.cpus}

    # Store Matplotlib cache inside the work dir to avoid $HOME permission issues
    export MPLCONFIGDIR=\$PWD/.mpl_cache

    plotCorrelation \
        -in  summary_counts.npz \
        --corMethod  spearman \
        --whatToPlot heatmap \
        --plotNumbers \
        --colorMap   coolwarm \
        --outFileCorMatrix  correlation_matrix.tsv \
        --plotFile         correlation_heatmap.png
    """
}

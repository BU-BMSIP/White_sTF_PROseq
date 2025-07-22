configfile: "config.yaml"

# Get sample names and paths from config
SAMPLES = config["samples"]
PROJECT_DIR = config["project_dir"]
REF_DIR = config["ref_dir"]

# Create output directories
rule all:
    input:
        # Quality control on trimmed files
        expand("results/fastqc_clean/{sample}_{read}_001_trimmed_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        
        # Quantification results (skip UMI extraction, go straight to quantification)
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLES),
        "results/featurecounts/gene_counts.txt",
        
        # Final analysis
        "results/deseq2/analysis_report.html"

# Quality control on trimmed files  
rule fastqc_clean:
    input:
        "rnaseq_data/trimmed/{sample}_{read}_001_trimmed.fastq.gz"
    output:
        html="results/fastqc_clean/{sample}_{read}_001_trimmed_fastqc.html",
        zip="results/fastqc_clean/{sample}_{read}_001_trimmed_fastqc.zip"
    threads: 2
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        "fastqc {input} -o results/fastqc_clean/ -t {threads}"

# Build Salmon index (if it doesn't exist)
rule build_salmon_index:
    input:
        fasta=config["genome_fasta"],
        gtf=config["gtf_file"]
    output:
        directory(config["salmon_index"])
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120
    shell:
        """
        salmon index \
            -t {input.fasta} \
            -i {output} \
            --gencode \
            -p {threads}
        """

# Salmon quantification (using trimmed files directly)
rule salmon_quant:
    input:
        r1="rnaseq_data/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        r2="rnaseq_data/trimmed/{sample}_R2_001_trimmed.fastq.gz",
        index=config["salmon_index"]
    output:
        quant="results/salmon/{sample}/quant.sf",
        lib="results/salmon/{sample}/lib_format_counts.json"
    params:
        outdir="results/salmon/{sample}"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240
    log:
        "logs/salmon/{sample}.log"
    shell:
        """
        salmon quant \
            -i {input.index} \
            -l A \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            --validateMappings \
            --gcBias \
            --seqBias \
            -p {threads} \
            &> {log}
        """

# Build STAR index (if it doesn't exist)
rule build_star_index:
    input:
        fasta=config["genome_fasta"],
        gtf=config["gtf_file"]
    output:
        directory(config["star_index"])
    threads: 16
    resources:
        mem_mb=32000,
        runtime=180
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang 100 \
             --runThreadN {threads}
        """

# STAR alignment (using trimmed files directly)
rule star_align:
    input:
        r1="rnaseq_data/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        r2="rnaseq_data/trimmed/{sample}_R2_001_trimmed.fastq.gz",
        index=config["star_index"]
    output:
        bam="results/star_alignment/{sample}_Aligned.sortedByCoord.out.bam",
        log_final="results/star_alignment/{sample}_Log.final.out"
    params:
        prefix="results/star_alignment/{sample}_"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=240
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard
        """

# FeatureCounts quantification
rule featurecounts:
    input:
        bams=expand("results/star_alignment/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        gtf=config["gtf_file"]
    output:
        counts="results/featurecounts/gene_counts.txt",
        summary="results/featurecounts/gene_counts.txt.summary"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120
    log:
        "logs/featurecounts.log"
    shell:
        """
        featureCounts \
            -T {threads} \
            -p \
            -B \
            -C \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} \
            &> {log}
        """

# DESeq2 analysis comparing both methods
rule deseq2_analysis:
    input:
        salmon_files=expand("results/salmon/{sample}/quant.sf", sample=SAMPLES),
        featurecounts="results/featurecounts/gene_counts.txt",
        design=config["experimental_design"]
    output:
        report="results/deseq2/analysis_report.html"
    resources:
        mem_mb=8000,
        runtime=60
    script:
        "scripts/deseq2_comparison.R"
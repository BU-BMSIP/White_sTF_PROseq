#!/bin/bash
# Spike-in alignment script for PRO-seq data

# Set up paths
PROJECT_DIR="/projectnb/khalil/nwhite42/ProSEQ_project"
DATA_DIR="${PROJECT_DIR}/data"
REF_DIR="${PROJECT_DIR}/refs/dm6"
OUTPUT_DIR="${PROJECT_DIR}/spike_in_counts"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check if bowtie2 is available
if ! command -v bowtie2 &> /dev/null; then
    echo "bowtie2 not found. Loading module..."
    module load bowtie2
fi

# Check if samtools is available
if ! command -v samtools &> /dev/null; then
    echo "samtools not found. Loading module..."
    module load samtools
fi

# Check if bowtie2 index exists
if [ ! -f "${REF_DIR}/drosophila_dm6.1.bt2" ]; then
    echo "Building bowtie2 index for dm6..."
    cd ${REF_DIR}
    # Assuming you have drosophila_dm6.fa or similar
    bowtie2-build drosophila_dm6.fa drosophila_dm6
    cd -
fi

echo "Starting spike-in alignments..."
echo "Sample,Total_reads,Spike_in_reads,Spike_in_percentage" > ${OUTPUT_DIR}/spike_in_summary.csv

# Process each fastq file
for fastq in ${DATA_DIR}/*.fastq.gz; do
    if [ -f "$fastq" ]; then
        base=$(basename $fastq .fastq.gz)
        echo "Processing $base..."
        
        # Count total reads in fastq
        echo "  Counting total reads..."
        total_reads=$(zcat $fastq | echo $((`wc -l`/4)))
        
        # Align to Drosophila
        echo "  Aligning to dm6 spike-in..."
        bowtie2 -x ${REF_DIR}/drosophila_dm6 \
                -U $fastq \
                -p 8 \
                --no-unal \
                -S ${OUTPUT_DIR}/${base}_dm6.sam \
                2> ${OUTPUT_DIR}/${base}_dm6_stats.txt
        
        # Count aligned reads
        spike_reads=$(samtools view -c -F 4 ${OUTPUT_DIR}/${base}_dm6.sam)
        spike_percent=$(echo "scale=2; $spike_reads * 100 / $total_reads" | bc)
        
        echo "  Total reads: $total_reads"
        echo "  Spike-in reads: $spike_reads ($spike_percent%)"
        
        # Add to summary
        echo "$base,$total_reads,$spike_reads,$spike_percent" >> ${OUTPUT_DIR}/spike_in_summary.csv
        
        # Remove SAM file to save space (optional)
        rm ${OUTPUT_DIR}/${base}_dm6.sam
    fi
done

echo "Done! Results saved to ${OUTPUT_DIR}/spike_in_summary.csv"

# Create R script for normalization
cat > ${OUTPUT_DIR}/apply_normalization.R << 'EOF'
# Read spike-in counts
spike_data <- read.csv("spike_in_summary.csv")
print(spike_data)

# Calculate normalization factors (higher spike-in = lower factor)
spike_data$norm_factor <- median(spike_data$Spike_in_reads) / spike_data$Spike_in_reads

# Save normalization factors
write.csv(spike_data, "spike_in_normalization_factors.csv", row.names = FALSE)

print("Normalization factors calculated. Use these to normalize your count matrix.")
print(spike_data[, c("Sample", "norm_factor")])
EOF

echo "Run the R script in ${OUTPUT_DIR} to calculate normalization factors"
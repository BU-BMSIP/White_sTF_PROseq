#!/usr/bin/env bash

echo "=== Converting Sample 19 BAM to BED and Generating Traces ==="

# Define paths
PROJECT_DIR="/projectnb/khalil/nwhite42/ProSEQ_project"
BAM_FILE="${PROJECT_DIR}/results/aligned/PR1570_S19_R1.sorted.bam"
BED_DIR="${PROJECT_DIR}/bed"
TRACES_DIR="${PROJECT_DIR}/traces"

# Create directories if they don't exist
mkdir -p "${BED_DIR}" "${TRACES_DIR}"

# Check if BAM file exists
if [ ! -f "${BAM_FILE}" ]; then
    echo "ERROR: BAM file not found: ${BAM_FILE}"
    echo "Available BAM files:"
    ls -la "${PROJECT_DIR}/results/aligned/"*S19*
    exit 1
fi

echo "Found BAM file: ${BAM_FILE}"
echo "File size: $(ls -lh ${BAM_FILE} | awk '{print $5}')"

# Load required modules
echo "Loading bedtools module..."
module load bedtools

# Convert BAM to BED
echo "Converting BAM to BED format..."
BED_OUTPUT="${BED_DIR}/PR1570_S19.bed"

bedtools bamtobed -i "${BAM_FILE}" > "${BED_OUTPUT}"

if [ -f "${BED_OUTPUT}" ]; then
    echo "✅ Successfully created BED file: ${BED_OUTPUT}"
    echo "Number of reads: $(wc -l < ${BED_OUTPUT})"
    echo "First few lines:"
    head -3 "${BED_OUTPUT}"
else
    echo "❌ Failed to create BED file"
    exit 1
fi

# Now run the trace generation for sample 19
echo ""
echo "=== Generating Traces and Metrics for Sample 19 ==="

python3 - "pSC146:" "8235" "8949" "10" <<'PY'
import sys, pandas as pd, numpy as np, pathlib
import warnings
warnings.filterwarnings('ignore')

contig = sys.argv[1]
tss, tes, bin_sz = map(int, sys.argv[2:5])

print(f"Processing for contig: '{contig}'")
print(f"TSS: {tss}, TES: {tes}, Bin size: {bin_sz}")

# Process the newly created PR1570_S19.bed file
bedfile = pathlib.Path('bed/PR1570_S19.bed')

if not bedfile.exists():
    print(f"ERROR: {bedfile} not found!")
    sys.exit(1)

print(f"Processing {bedfile.name}...")

# Read BED file
try:
    df = pd.read_csv(bedfile, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'name', 'score', 'strand'])
except Exception as e:
    print(f"Error reading {bedfile}: {e}")
    sys.exit(1)

print(f"  Total reads: {len(df)}")
print(f"  Available chromosomes: {sorted(df['chr'].unique()[:10])}")  # Show first 10 chroms

# Filter for target contig (exact match)
df_contig = df[df['chr'] == contig].copy()
print(f"  Reads on '{contig}': {len(df_contig)}")

if len(df_contig) == 0:
    print(f"  No reads found for contig '{contig}'")
    print(f"  Available contigs (first 20): {sorted(df['chr'].unique())[:20]}")
    
    # Let's try some common alternatives
    alternatives = ['pSC146', 'chr1', 'Chr1', '1']
    for alt in alternatives:
        alt_count = len(df[df['chr'] == alt])
        if alt_count > 0:
            print(f"  Found {alt_count} reads on '{alt}' - you might want to use this instead")
    sys.exit(1)

# Extract 3' end positions (single nucleotide where Pol II stopped)
# For PRO-seq: + strand reads use start position, - strand reads use end-1
df_contig['three_prime'] = np.where(df_contig['strand'] == '+', 
                                   df_contig['start'], 
                                   df_contig['end'] - 1)

# Count reads at each position
counts = df_contig['three_prime'].value_counts().sort_index()

if len(counts) == 0:
    print(f"  No valid 3' positions found")
    sys.exit(1)

print(f"  Position range: {counts.index.min()} - {counts.index.max()}")

# Create binned counts
min_pos = counts.index.min()
max_pos = counts.index.max()
bin_edges = range(min_pos - (min_pos % bin_sz), max_pos + bin_sz, bin_sz)

binned_counts = []
for i in range(len(bin_edges) - 1):
    bin_start = bin_edges[i]
    bin_end = bin_edges[i + 1]
    bin_count = counts.loc[bin_start:bin_end-1].sum()
    if bin_count > 0:
        binned_counts.append((bin_start, bin_count))

# Save binned data
out_base = pathlib.Path('traces') / bedfile.stem

if binned_counts:
    binned_df = pd.DataFrame(binned_counts, columns=['bin_start', 'count'])
    binned_df.to_csv(f'{out_base}.binned.tsv', sep='\t', header=False, index=False)
    print(f"  Saved binned data: {len(binned_counts)} bins")
else:
    print(f"  No binned data to save")

# Calculate pausing index and metrics
try:
    # Pausing index: ratio of reads near TSS vs downstream
    tss_region = counts.loc[max(0, tss-50):tss+50].sum()
    body_region = counts.loc[tss+50:min(tss+250, max_pos)].sum()
    
    if body_region > 0:
        pausing_index = tss_region / body_region
    else:
        pausing_index = float('inf')
    
    # Regional analysis
    regions = [
        ('promoter', max(0, tss-50), tss+50),
        ('gene_body', tss+50, tes-150),
        ('pre_TES', tes-150, tes),
        ('post_TES', tes, tes+150),
        ('downstream', tes+150, min(tes+500, max_pos))
    ]
    
    metrics = [f'pausing_index\t{pausing_index:.4f}']
    
    for region_name, start, end in regions:
        if start < end and start <= max_pos:
            region_sum = counts.loc[start:min(end, max_pos)].sum()
            region_length = min(end, max_pos) - start + 1
            region_density = region_sum / region_length if region_length > 0 else 0
            metrics.append(f'{region_name}_density\t{region_density:.4f}')
            metrics.append(f'{region_name}_total\t{region_sum}')
    
    # Save metrics
    with open(f'{out_base}.metrics.txt', 'w') as f:
        f.write('\n'.join(metrics) + '\n')
    
    print(f"  Pausing index: {pausing_index:.4f}")
    
except Exception as e:
    print(f"  Error calculating metrics: {e}")

print("Processing complete for PR1570_S19!")
PY

echo ""
echo "=== Sample 19 Processing Complete ==="
echo "Files created:"
echo "  ${BED_OUTPUT}"
echo "  ${TRACES_DIR}/PR1570_S19.binned.tsv"
echo "  ${TRACES_DIR}/PR1570_S19.metrics.txt"
echo ""
echo "Verifying results:"
ls -la "${BED_OUTPUT}" "${TRACES_DIR}/PR1570_S19."*

echo ""
echo "=== All samples now processed ==="
echo "BED files available:"
ls -1 bed/*.bed | sed 's/.*\///' | sed 's/\.bed$//'

echo ""
echo "Trace files available:"
ls -1 traces/*.binned.tsv 2>/dev/null | sed 's/.*\///' | sed 's/\.binned\.tsv$//' || echo "No trace files found yet"
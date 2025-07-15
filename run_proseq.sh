#!/usr/bin/env bash
set -euo pipefail

#################################  USAGE  #################################
#   bash run_proseq.sh -s samples_R1_only.csv \
#                      -f refs/merged_reference.fa \
#                      -c 'pSC146' \
#                      -t 8235 -e 8949
############################################################################

# defaults ---------------------------------------------------------
CSV=''                    # sample sheet
REF_FA=''                 # merged reference fasta
CONTIG=''                 # plasmid header (without colon)
TSS=0 TES=0               # coordinates (0-based)
BIN=10
THREADS=8                 # change if you want bowtie2/samtools faster
ADAPTER='TGGAATTCTCGG'
UMI_PATTERN='NNNNNN'

while getopts "s:f:c:t:e:b:p:" opt; do
  case $opt in
    s) CSV=$OPTARG ;;
    f) REF_FA=$OPTARG ;;
    c) CONTIG=$OPTARG ;;
    t) TSS=$OPTARG ;;
    e) TES=$OPTARG ;;
    b) BIN=$OPTARG ;;
    p) THREADS=$OPTARG ;;
  esac
done

if [[ -z $CSV || -z $REF_FA || -z $CONTIG ]]; then
  echo "Usage: bash run_proseq.sh -s samples.csv -f ref.fa -c plasmidName -t TSS -e TES"
  echo "Example: bash run_proseq.sh -s samples.csv -f ref.fa -c pSC146 -t 8235 -e 8949"
  exit 1
fi

# Remove trailing colon if present
CONTIG=${CONTIG%:}

echo "Starting PRO-seq analysis..."
echo "Reference: $REF_FA"
echo "Contig: $CONTIG"
echo "TSS: $TSS, TES: $TES"

mkdir -p trimmed umi bam dedup bed traces

# Check if CSV file exists and has correct format
if [[ ! -f "$CSV" ]]; then
  echo "Error: Sample sheet $CSV not found!"
  exit 1
fi

echo "Sample sheet preview:"
head -3 "$CSV"

# 1. Adapter trim + UMI extract -----------------------------------
echo "=== Step 1: Adapter trimming and UMI extraction ==="
tail -n +2 "$CSV" | while IFS=',' read -r ID FQ ; do
  # Remove any quotes and whitespace
  ID=$(echo "$ID" | tr -d '"' | tr -d ' ')
  FQ=$(echo "$FQ" | tr -d '"' | tr -d ' ')
  
  echo "Processing sample: $ID"
  echo "  Input file: $FQ"
  
  if [[ ! -f "$FQ" ]]; then
    echo "  WARNING: File $FQ not found, skipping..."
    continue
  fi
  
  # Adapter trimming
  cutadapt -a "$ADAPTER" -m 26 -q 10 -o trimmed/${ID}.fq.gz "$FQ" \
    --report=minimal
  
  # UMI extraction
  umi_tools extract --stdin trimmed/${ID}.fq.gz  \
                    --stdout umi/${ID}.umi.fq.gz \
                    --bc-pattern="$UMI_PATTERN" \
                    --log2stderr
done

# 2. Build bowtie2 index once -------------------------------------
echo "=== Step 2: Building Bowtie2 index ==="
if [[ ! -e ${REF_FA}.1.bt2 ]]; then
  echo "Building Bowtie2 index for $REF_FA..."
  bowtie2-build "$REF_FA" "${REF_FA%.fa}"
else
  echo "Bowtie2 index already exists, skipping..."
fi

# 3. Align, sort, index -------------------------------------------
echo "=== Step 3: Alignment and sorting ==="
for FQ in umi/*.umi.fq.gz ; do
  if [[ ! -f "$FQ" ]]; then
    echo "No UMI files found, check previous steps"
    continue
  fi
  
  ID=$(basename "$FQ" .umi.fq.gz)
  echo "Aligning $ID..."
  
  bowtie2 -x "${REF_FA%.fa}" -U "$FQ" --very-sensitive -p "$THREADS" \
    --no-unal 2>/dev/null | \
    samtools sort -@ "$THREADS" -o bam/${ID}.bam
  
  samtools index bam/${ID}.bam
  
  # Quick alignment stats
  echo "  Alignment stats for $ID:"
  samtools flagstat bam/${ID}.bam | head -1
done

# 4. Optional UMI dedup; comment out if not wanted ----------------
echo "=== Step 4: UMI deduplication ==="
for B in bam/*.bam ; do
  if [[ ! -f "$B" ]]; then continue; fi
  
  ID=$(basename "$B" .bam)
  echo "Deduplicating $ID..."
  
  umi_tools dedup -I "$B" -S dedup/${ID}.bam --log2stderr
done

# If you want to skip dedup, comment out the above and uncomment:
# ln -sf ../bam/*.bam dedup/

# 5. BAM -> BED ----------------------------------------------------
echo "=== Step 5: Converting BAM to BED ==="
for B in dedup/*.bam ; do
  if [[ ! -f "$B" ]]; then continue; fi
  
  ID=$(basename "$B" .bam)
  echo "Converting $ID to BED..."
  
  bedtools bamtobed -i "$B" > bed/${ID}.bed
  
  # Show some stats
  echo "  Total reads in BED: $(wc -l < bed/${ID}.bed)"
  echo "  Reads on $CONTIG: $(grep -c "^${CONTIG}" bed/${ID}.bed || echo 0)"
done

# 6. 3'-end counting / pausing index ------------------------------
echo "=== Step 6: Generating traces and metrics ==="
python3 - "$CONTIG" "$TSS" "$TES" "$BIN" <<'PY'
import sys, pandas as pd, numpy as np, pathlib
import warnings
warnings.filterwarnings('ignore')

contig = sys.argv[1]
tss, tes, bin_sz = map(int, sys.argv[2:5])

print(f"Processing for contig: {contig}")
print(f"TSS: {tss}, TES: {tes}, Bin size: {bin_sz}")

bed_files = list(pathlib.Path('bed').glob('*.bed'))
if not bed_files:
    print("No BED files found!")
    sys.exit(1)

for bedfile in bed_files:
    print(f"Processing {bedfile.name}...")
    
    # Read BED file
    try:
        df = pd.read_csv(bedfile, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'name', 'score', 'strand'])
    except Exception as e:
        print(f"Error reading {bedfile}: {e}")
        continue
    
    print(f"  Total reads: {len(df)}")
    
    # Filter for target contig
    df_contig = df[df['chr'] == contig].copy()
    print(f"  Reads on {contig}: {len(df_contig)}")
    
    if len(df_contig) == 0:
        print(f"  No reads found for contig {contig}")
        continue
    
    # Extract 3' end positions (single nucleotide where Pol II stopped)
    # For PRO-seq: + strand reads use start position, - strand reads use end-1
    df_contig['three_prime'] = np.where(df_contig['strand'] == '+', 
                                       df_contig['start'], 
                                       df_contig['end'] - 1)
    
    # Count reads at each position
    counts = df_contig['three_prime'].value_counts().sort_index()
    
    if len(counts) == 0:
        print(f"  No valid 3' positions found")
        continue
    
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

print("Processing complete!")
PY

echo ""
echo "=== Analysis Complete ==="
echo "Results saved in:"
echo "  traces/*.binned.tsv   - Binned read counts"
echo "  traces/*.metrics.txt  - Pausing indices and regional metrics"
echo ""
echo "Quick plot command:"
echo "python3 -c \"import pandas as pd, matplotlib.pyplot as plt; df=pd.read_csv('traces/$(ls traces/*.binned.tsv | head -1 | xargs basename)', sep='\t', header=None, names=['pos','count']); plt.figure(figsize=(12,4)); plt.plot(df.pos, df.count); plt.axvline($TSS, color='red', linestyle='--', label='TSS'); plt.axvline($TES, color='blue', linestyle='--', label='TES'); plt.xlabel('Position'); plt.ylabel('Read Count'); plt.legend(); plt.title('PRO-seq Trace'); plt.tight_layout(); plt.show()\""
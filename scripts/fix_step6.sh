#!/usr/bin/env bash

# Quick fix - rerun only Step 6 with correct contig name
echo "=== Step 6 (FIXED): Generating traces and metrics ==="

python3 - "pSC146:" "8235" "8949" "10" <<'PY'
import sys, pandas as pd, numpy as np, pathlib
import warnings
warnings.filterwarnings('ignore')

contig = sys.argv[1]
tss, tes, bin_sz = map(int, sys.argv[2:5])

print(f"Processing for contig: '{contig}'")
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
    
    # Filter for target contig (exact match)
    df_contig = df[df['chr'] == contig].copy()
    print(f"  Reads on '{contig}': {len(df_contig)}")
    
    if len(df_contig) == 0:
        print(f"  No reads found for contig '{contig}'")
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
echo "=== Fixed Analysis Complete ==="
echo "Results saved in:"
echo "  traces/*.binned.tsv   - Binned read counts"
echo "  traces/*.metrics.txt  - Pausing indices and regional metrics"
echo ""
echo "Check results:"
ls -la traces/
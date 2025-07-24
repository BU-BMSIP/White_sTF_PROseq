#!/usr/bin/env bash

echo "=== CORRECTED PRO-seq Analysis: Minus Strand Only ==="

python3 - "pSC146:" "8235" "8949" "10" <<'PY'
import sys, pandas as pd, numpy as np, pathlib
import warnings
warnings.filterwarnings('ignore')

contig = sys.argv[1]
tss, tes, bin_sz = map(int, sys.argv[2:5])

print(f"Processing for contig: '{contig}'")
print(f"TSS: {tss}, TES: {tes}, Bin size: {bin_sz}")
print("*** USING CORRECTED PRO-seq METHODOLOGY ***")
print("*** R1 reads on MINUS strand only (3' end positions) ***")

bed_files = list(pathlib.Path('bed').glob('*.bed'))
if not bed_files:
    print("No BED files found!")
    sys.exit(1)

for bedfile in bed_files:
    print(f"\nProcessing {bedfile.name}...")
    
    # Read BED file
    try:
        df = pd.read_csv(bedfile, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'name', 'score', 'strand'])
    except Exception as e:
        print(f"Error reading {bedfile}: {e}")
        continue
    
    print(f"  Total reads: {len(df)}")
    print(f"  Plus strand reads: {len(df[df['strand'] == '+'])}")
    print(f"  Minus strand reads: {len(df[df['strand'] == '-'])}")
    
    # Filter for target contig (exact match)
    df_contig = df[df['chr'] == contig].copy()
    print(f"  Reads on '{contig}': {len(df_contig)}")
    
    if len(df_contig) == 0:
        print(f"  No reads found for contig '{contig}'")
        # Show available contigs for debugging
        available_contigs = sorted(df['chr'].unique())[:10]
        print(f"  Available contigs (first 10): {available_contigs}")
        continue
    
    # *** CRITICAL: PRO-seq uses ONLY minus strand reads ***
    # These represent R1 reads where the original nascent RNA was on the plus strand
    df_minus = df_contig[df_contig['strand'] == '-'].copy()
    print(f"  Minus strand reads on '{contig}': {len(df_minus)}")
    
    if len(df_minus) == 0:
        print(f"  No minus strand reads found for contig '{contig}'")
        continue
    
    # For minus strand reads: the 3' end position is the START coordinate
    # This represents where Pol II was when it incorporated the biotin-labeled nucleotide
    df_minus['pol2_position'] = df_minus['start']
    
    # Count reads at each Pol II position
    counts = df_minus['pol2_position'].value_counts().sort_index()
    
    if len(counts) == 0:
        print(f"  No valid Pol II positions found")
        continue
    
    print(f"  Position range: {counts.index.min()} - {counts.index.max()}")
    print(f"  Unique positions: {len(counts)}")
    print(f"  Total signal: {counts.sum()}")
    
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
    
    # Save corrected binned data
    out_base = pathlib.Path('traces') / f"{bedfile.stem}_corrected"
    
    if binned_counts:
        binned_df = pd.DataFrame(binned_counts, columns=['bin_start', 'count'])
        binned_df.to_csv(f'{out_base}.binned.tsv', sep='\t', header=False, index=False)
        print(f"  Saved corrected binned data: {len(binned_counts)} bins")
    
    # Calculate pausing index and metrics
    try:
        # Pausing index: ratio of reads near TSS vs downstream gene body
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
        metrics.append(f'total_minus_strand_reads\t{len(df_minus)}')
        metrics.append(f'total_all_reads\t{len(df_contig)}')
        metrics.append(f'minus_strand_fraction\t{len(df_minus)/len(df_contig):.4f}')
        
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
        
        print(f"  Pausing index (corrected): {pausing_index:.4f}")
        print(f"  Minus strand fraction: {len(df_minus)/len(df_contig):.3f}")
        
    except Exception as e:
        print(f"  Error calculating metrics: {e}")

print("\n*** CORRECTED PRO-seq Analysis Complete! ***")
print("Results saved with '_corrected' suffix:")
print("  traces/*_corrected.binned.tsv   - Corrected binned read counts (minus strand only)")
print("  traces/*_corrected.metrics.txt  - Corrected pausing indices and metrics")
PY

echo ""
echo "=== Comparison with Original Analysis ==="
echo "Compare file sizes - corrected should be smaller (minus strand only):"
echo "Original files:"
ls -la traces/*.binned.tsv | grep -v corrected
echo ""
echo "Corrected files:"
ls -la traces/*_corrected.binned.tsv 2>/dev/null || echo "No corrected files yet"
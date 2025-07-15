#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np

# Plot comparison: Original vs Corrected traces
fig, axes = plt.subplots(4, 2, figsize=(20, 16))

# Get files
original_files = sorted(glob.glob('traces/*.binned.tsv'))
corrected_files = sorted(glob.glob('traces/*_corrected.binned.tsv'))

TSS, TES = 8235, 8949

print(f"Found {len(original_files)} original files")
print(f"Found {len(corrected_files)} corrected files")

# Filter original files to exclude corrected ones
original_files = [f for f in original_files if '_corrected' not in f]

print(f"Comparing {len(original_files)} original vs {len(corrected_files)} corrected")

for i in range(min(4, len(original_files))):
    # Original trace
    if i < len(original_files):
        orig_file = original_files[i]
        sample_id = orig_file.split('/')[-1].replace('.binned.tsv', '')
        
        df_orig = pd.read_csv(orig_file, sep='\t', header=None, names=['pos', 'count'])
        
        ax_orig = axes[i, 0]
        ax_orig.plot(df_orig.pos, df_orig['count'], linewidth=1.5, color='steelblue', label=f'{sample_id} (Original)')
        ax_orig.axvline(TSS, color='red', linestyle='--', alpha=0.7, label='TSS')
        ax_orig.axvline(TES, color='blue', linestyle='--', alpha=0.7, label='TES')
        ax_orig.set_title(f'{sample_id} - Original (Both Strands)')
        ax_orig.set_xlabel('Position (bp)')
        ax_orig.set_ylabel('Read Count')
        ax_orig.legend()
        ax_orig.grid(True, alpha=0.3)
        ax_orig.set_xlim(7000, 11000)
        
        # Add stats
        total_reads_orig = df_orig['count'].sum()
        max_signal_orig = df_orig['count'].max()
        ax_orig.text(0.02, 0.98, f'Total: {total_reads_orig:,}\nMax: {max_signal_orig:,}', 
                    transform=ax_orig.transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Corrected trace
    if i < len(corrected_files):
        corr_file = corrected_files[i]
        sample_id_corr = corr_file.split('/')[-1].replace('_corrected.binned.tsv', '')
        
        df_corr = pd.read_csv(corr_file, sep='\t', header=None, names=['pos', 'count'])
        
        ax_corr = axes[i, 1]
        ax_corr.plot(df_corr.pos, df_corr['count'], linewidth=1.5, color='darkgreen', label=f'{sample_id_corr} (Corrected)')
        ax_corr.axvline(TSS, color='red', linestyle='--', alpha=0.7, label='TSS')
        ax_corr.axvline(TES, color='blue', linestyle='--', alpha=0.7, label='TES')
        ax_corr.set_title(f'{sample_id_corr} - Corrected (Minus Strand Only)')
        ax_corr.set_xlabel('Position (bp)')
        ax_corr.set_ylabel('Read Count')
        ax_corr.legend()
        ax_corr.grid(True, alpha=0.3)
        ax_corr.set_xlim(7000, 11000)
        
        # Add stats
        total_reads_corr = df_corr['count'].sum()
        max_signal_corr = df_corr['count'].max()
        ax_corr.text(0.02, 0.98, f'Total: {total_reads_corr:,}\nMax: {max_signal_corr:,}', 
                    transform=ax_corr.transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

plt.suptitle('PRO-seq Traces: Original vs Corrected Methodology', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('original_vs_corrected_traces.png', dpi=300, bbox_inches='tight')
plt.show()

print("Comparison plot saved as: original_vs_corrected_traces.png")

# Generate detailed comparison statistics
print("\n" + "="*60)
print("DETAILED COMPARISON STATISTICS")
print("="*60)

comparison_data = []

for i in range(min(len(original_files), len(corrected_files))):
    orig_file = original_files[i]
    corr_file = corrected_files[i]
    
    sample_id = orig_file.split('/')[-1].replace('.binned.tsv', '')
    
    df_orig = pd.read_csv(orig_file, sep='\t', header=None, names=['pos', 'count'])
    df_corr = pd.read_csv(corr_file, sep='\t', header=None, names=['pos', 'count'])
    
    # Calculate statistics
    orig_total = df_orig['count'].sum()
    corr_total = df_corr['count'].sum()
    
    orig_max = df_orig['count'].max()
    corr_max = df_corr['count'].max()
    
    # TSS region analysis
    tss_region_orig = df_orig[(df_orig.pos >= TSS-50) & (df_orig.pos <= TSS+50)]['count'].sum()
    tss_region_corr = df_corr[(df_corr.pos >= TSS-50) & (df_corr.pos <= TSS+50)]['count'].sum()
    
    # Signal reduction
    total_reduction = (orig_total - corr_total) / orig_total * 100
    tss_reduction = (tss_region_orig - tss_region_corr) / tss_region_orig * 100 if tss_region_orig > 0 else 0
    
    print(f"\n{sample_id}:")
    print(f"  Total reads:     {orig_total:,} → {corr_total:,} ({total_reduction:.1f}% reduction)")
    print(f"  Max signal:      {orig_max:,} → {corr_max:,}")
    print(f"  TSS region:      {tss_region_orig:,} → {tss_region_corr:,} ({tss_reduction:.1f}% reduction)")
    print(f"  Signal ratio:    {corr_total/orig_total:.3f}")
    
    comparison_data.append({
        'Sample': sample_id,
        'Original_Total': orig_total,
        'Corrected_Total': corr_total,
        'Reduction_Percent': total_reduction,
        'Signal_Ratio': corr_total/orig_total,
        'Original_Max': orig_max,
        'Corrected_Max': corr_max
    })

# Create summary DataFrame
summary_df = pd.DataFrame(comparison_data)
summary_df.to_csv('original_vs_corrected_comparison.csv', index=False)

print(f"\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Average signal reduction: {summary_df['Reduction_Percent'].mean():.1f}%")
print(f"Average signal ratio: {summary_df['Signal_Ratio'].mean():.3f}")
print(f"Expected reduction: ~50% (if strands are balanced)")
print(f"Summary saved to: original_vs_corrected_comparison.csv")

# Check if the reduction is as expected
expected_reduction = 50
actual_reduction = summary_df['Reduction_Percent'].mean()
if abs(actual_reduction - expected_reduction) < 10:
    print(f"✅ Signal reduction ({actual_reduction:.1f}%) is close to expected (~{expected_reduction}%)")
else:
    print(f"⚠️  Signal reduction ({actual_reduction:.1f}%) differs from expected (~{expected_reduction}%)")
    print("   This could indicate strand bias in your data.")
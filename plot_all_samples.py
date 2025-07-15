#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import glob

# Simple diagnostic plot for all 4 samples
fig, axes = plt.subplots(4, 1, figsize=(15, 16))

files = glob.glob('traces/*.binned.tsv')
files.sort()

TSS, TES = 8235, 8949

print(f"Found {len(files)} trace files:")
for f in files:
    print(f"  {f}")
print()

for i, file in enumerate(files):
    sample_id = file.split('/')[-1].replace('.binned.tsv', '')
    
    # Read data
    df = pd.read_csv(file, sep='\t', header=None, names=['pos', 'count'])
    
    print(f"Sample {sample_id}:")
    print(f"  Total positions: {len(df)}")
    print(f"  Position range: {df.pos.min()} - {df.pos.max()}")
    print(f"  Count range: {df['count'].min()} - {df['count'].max()}")
    print(f"  Total reads: {df['count'].sum()}")
    
    # Find peak signal
    max_idx = df['count'].idxmax()
    print(f"  Max signal: {df.loc[max_idx, 'count']} at position {df.loc[max_idx, 'pos']}")
    
    # Check signal around TSS and TES
    tss_region = df[(df.pos >= TSS-100) & (df.pos <= TSS+100)]
    tes_region = df[(df.pos >= TES-100) & (df.pos <= TES+100)]
    
    if len(tss_region) > 0:
        print(f"  TSS region (±100bp): max = {tss_region['count'].max()}")
    if len(tes_region) > 0:
        print(f"  TES region (±100bp): max = {tes_region['count'].max()}")
    print()
    
    # Plot
    ax = axes[i]
    ax.plot(df.pos, df['count'], linewidth=1.5, label=sample_id)
    ax.axvline(TSS, color='red', linestyle='--', alpha=0.7, label='TSS')
    ax.axvline(TES, color='blue', linestyle='--', alpha=0.7, label='TES')
    ax.set_title(f'{sample_id} - Raw PRO-seq Signal')
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Read Count')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Focus on gene region
    ax.set_xlim(7000, 11000)

plt.tight_layout()
plt.savefig('diagnostic_traces_all_samples.png', dpi=300, bbox_inches='tight')
plt.show()

print("Diagnostic plot saved as: diagnostic_traces_all_samples.png")

# Summary statistics
print("\n=== SUMMARY STATISTICS ===")
summary_data = []
for file in files:
    sample_id = file.split('/')[-1].replace('.binned.tsv', '')
    df = pd.read_csv(file, sep='\t', header=None, names=['pos', 'count'])
    
    # Calculate key metrics
    total_reads = df['count'].sum()
    max_signal = df['count'].max()
    max_pos = df.loc[df['count'].idxmax(), 'pos']
    
    # TSS and TES region signals
    tss_region = df[(df.pos >= TSS-50) & (df.pos <= TSS+50)]
    tes_region = df[(df.pos >= TES-50) & (df.pos <= TES+50)]
    gene_body = df[(df.pos >= TSS+50) & (df.pos <= TES-50)]
    
    tss_signal = tss_region['count'].sum() if len(tss_region) > 0 else 0
    tes_signal = tes_region['count'].sum() if len(tes_region) > 0 else 0
    body_signal = gene_body['count'].sum() if len(gene_body) > 0 else 0
    
    # Pausing index (TSS/body ratio)
    pausing_index = tss_signal / body_signal if body_signal > 0 else float('inf')
    
    summary_data.append({
        'Sample': sample_id,
        'Total_Reads': total_reads,
        'Max_Signal': max_signal,
        'Max_Position': max_pos,
        'TSS_Signal': tss_signal,
        'TES_Signal': tes_signal,
        'Body_Signal': body_signal,
        'Pausing_Index': pausing_index
    })

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))

# Save summary
summary_df.to_csv('sample_summary_statistics.csv', index=False)
print(f"\nSummary statistics saved to: sample_summary_statistics.csv")
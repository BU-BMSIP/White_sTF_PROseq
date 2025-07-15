#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob

# Set up the plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10))

# Color scheme
colors = {
    'PR1567_S16': '#1f77b4',  # Blue for 656-BR1
    'PR1568_S17': '#ff7f0e',  # Orange for 917-BR1  
    'PR1569_S18': '#2ca02c',  # Green for 656-BR2
    'PR1570_S19': '#d62728'   # Red for 917-BR2 (if exists)
}

# Sample info
sample_info = {
    'PR1567_S16': '656-BR1 (Synthetic TF)',
    'PR1568_S17': '917-BR1 (Control)', 
    'PR1569_S18': '656-BR2 (Synthetic TF)',
    'PR1570_S19': '917-BR2 (Control)'
}

# TSS and TES positions
TSS = 8235
TES = 8949

# Plot 1: Individual traces
print("=== PRO-seq Trace Analysis ===")
print(f"TSS: {TSS}, TES: {TES}")
print()

binned_files = glob.glob('traces/*.binned.tsv')
binned_files.sort()

for file in binned_files:
    sample_id = file.split('/')[-1].replace('.binned.tsv', '')
    
    print(f"Processing {sample_id}...")
    
    # Read the data with debugging
    try:
        df = pd.read_csv(file, sep='\t', header=None)
        print(f"  Data shape: {df.shape}")
        print(f"  First few rows:")
        print(df.head())
        
        # Handle different possible formats
        if df.shape[1] == 2:
            df.columns = ['pos', 'count']
        else:
            print(f"  Unexpected format, skipping {sample_id}")
            continue
            
        # Basic validation
        if len(df) == 0:
            print(f"  No data in {sample_id}")
            continue
            
        # Plot individual trace
        ax1.plot(df.pos, df.count, 
                 color=colors.get(sample_id, 'gray'), 
                 label=sample_info.get(sample_id, sample_id),
                 linewidth=2, alpha=0.8)
    except Exception as e:
        print(f"  Error processing {sample_id}: {e}")
        continue
    
    print(f"Sample: {sample_id}")
    print(f"  Position range: {df.pos.min()} - {df.pos.max()}")
    print(f"  Total reads: {df['count'].sum()}")
    print(f"  Max signal: {df['count'].max()} at position {df.loc[df['count'].idxmax(), 'pos']}")
    
    # Read metrics
    metrics_file = file.replace('.binned.tsv', '.metrics.txt')
    try:
        with open(metrics_file, 'r') as f:
            metrics = dict(line.strip().split('\t') for line in f)
        print(f"  Pausing index: {float(metrics['pausing_index']):.3f}")
    except:
        print(f"  Could not read metrics from {metrics_file}")
    print()

# Formatting for plot 1
ax1.axvline(TSS, color='red', linestyle='--', alpha=0.7, label='TSS')
ax1.axvline(TES, color='blue', linestyle='--', alpha=0.7, label='TES')
ax1.set_xlabel('Plasmid Position (bp)')
ax1.set_ylabel('PRO-seq Signal (Read Count)')
ax1.set_title('PRO-seq Traces - Individual Samples')
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(7000, 11000)  # Focus on the gene region

# Plot 2: Average by condition (656 vs 917)
condition_data = {'656': [], '917': []}

for file in binned_files:
    sample_id = file.split('/')[-1].replace('.binned.tsv', '')
    df = pd.read_csv(file, sep='\t', header=None, names=['pos', 'count'])
    
    if '656' in sample_id:
        condition_data['656'].append(df)
    elif '917' in sample_id:
        condition_data['917'].append(df)

# Calculate averages
for condition, data_list in condition_data.items():
    if not data_list:
        continue
        
    # Merge all dataframes for this condition
    all_positions = set()
    for df in data_list:
        all_positions.update(df.pos.values)
    
    all_positions = sorted(all_positions)
    
    # Create average trace
    avg_counts = []
    for pos in all_positions:
        counts_at_pos = []
        for df in data_list:
            count = df[df.pos == pos]['count'].values
            counts_at_pos.append(count[0] if len(count) > 0 else 0)
        avg_counts.append(np.mean(counts_at_pos))
    
    # Plot average
    color = '#1f77b4' if condition == '656' else '#ff7f0e'
    label = f'{condition} (Synthetic TF, n={len(data_list)})' if condition == '656' else f'{condition} (Control, n={len(data_list)})'
    
    ax2.plot(all_positions, avg_counts, color=color, label=label, linewidth=3)
    
    # Add confidence interval if multiple replicates
    if len(data_list) > 1:
        std_counts = []
        for pos in all_positions:
            counts_at_pos = []
            for df in data_list:
                count = df[df.pos == pos]['count'].values
                counts_at_pos.append(count[0] if len(count) > 0 else 0)
            std_counts.append(np.std(counts_at_pos))
        
        ax2.fill_between(all_positions, 
                        np.array(avg_counts) - np.array(std_counts),
                        np.array(avg_counts) + np.array(std_counts),
                        color=color, alpha=0.2)

# Formatting for plot 2
ax2.axvline(TSS, color='red', linestyle='--', alpha=0.7, label='TSS')
ax2.axvline(TES, color='blue', linestyle='--', alpha=0.7, label='TES')
ax2.set_xlabel('Plasmid Position (bp)')
ax2.set_ylabel('PRO-seq Signal (Read Count)')
ax2.set_title('PRO-seq Traces - Average by Condition')
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xlim(7000, 11000)

plt.tight_layout()
plt.savefig('proseq_traces_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

print("=== Analysis Summary ===")
print("Plot saved as: proseq_traces_comparison.png")
print()
print("Key observations to look for:")
print("1. Pausing at TSS (red line) - higher values = more pausing")
print("2. Signal in gene body (between red and blue lines)")  
print("3. Termination pattern after TES (blue line)")
print("4. Overall signal differences between 656 vs 917 conditions")
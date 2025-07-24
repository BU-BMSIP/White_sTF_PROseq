#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
from scipy import stats
import seaborn as sns

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'legend.fontsize': 11,
    'figure.titlesize': 16
})

def load_and_process_traces(file_pattern, condition_name):
    """Load trace files and return combined dataframe"""
    files = sorted(glob.glob(file_pattern))
    print(f"Found {len(files)} files for {condition_name}")
    
    all_data = []
    for i, file in enumerate(files):
        df = pd.read_csv(file, sep='\t', header=None, names=['pos', 'count'])
        df['sample'] = f"{condition_name}_Rep{i+1}"
        df['condition'] = condition_name
        all_data.append(df)
    
    return pd.concat(all_data, ignore_index=True) if all_data else pd.DataFrame()

def create_binned_average(data, bin_size=10, position_range=None):
    """Create binned averages with confidence intervals"""
    if data.empty:
        return pd.DataFrame()
        
    if position_range:
        data = data[(data['pos'] >= position_range[0]) & (data['pos'] <= position_range[1])]
    
    # Create bins
    min_pos = data['pos'].min()
    max_pos = data['pos'].max()
    
    # Align bins to bin_size
    bin_start = (min_pos // bin_size) * bin_size
    bin_end = ((max_pos // bin_size) + 1) * bin_size
    
    bins = range(int(bin_start), int(bin_end), bin_size)
    
    # Assign bins to positions
    data['bin'] = pd.cut(data['pos'], bins=bins, labels=bins[:-1], include_lowest=True)
    data['bin'] = data['bin'].astype(float)
    
    # Calculate statistics for each bin
    bin_stats = data.groupby(['bin', 'sample'])['count'].sum().reset_index()
    
    # Calculate mean and confidence intervals
    summary = bin_stats.groupby('bin')['count'].agg([
        'mean', 'std', 'count', 'sem'
    ]).reset_index()
    
    # Calculate 95% confidence intervals
    summary['ci_lower'] = summary['mean'] - 1.96 * summary['sem']
    summary['ci_upper'] = summary['mean'] + 1.96 * summary['sem']
    
    # Make sure CI doesn't go below 0
    summary['ci_lower'] = np.maximum(summary['ci_lower'], 0)
    
    return summary

# Parameters
TSS, TES = 8235, 8949
PLOT_RANGE = (TSS - 500, TES + 500)  # 500bp upstream/downstream
BIN_SIZE = 10

print("=== Creating Multi-Condition PRO-seq Traces ===")

# Define your experimental conditions
conditions = {
    'p65': 'traces/PR1567_S16_corrected.binned.tsv traces/PR1568_S18_corrected.binned.tsv',
    '2xTIMs-p65': 'traces/PR1569_S17_corrected.binned.tsv traces/PR1570_S19_corrected.binned.tsv'
}

# Alternative: treat all as replicates of one condition
# conditions = {
#     'PRO-seq Signal': 'traces/*_corrected.binned.tsv'
# }

# Load data for each condition
condition_data = {}
condition_summaries = {}

for condition_name, file_pattern in conditions.items():
    print(f"\nProcessing {condition_name}...")
    
    # Handle space-separated file lists vs glob patterns
    if ' ' in file_pattern:
        files = file_pattern.split()
        all_data = []
        for i, file in enumerate(files):
            if glob.glob(file):  # Check if file exists
                df = pd.read_csv(file, sep='\t', header=None, names=['pos', 'count'])
                df['sample'] = f"{condition_name}_Rep{i+1}"
                df['condition'] = condition_name
                all_data.append(df)
        data = pd.concat(all_data, ignore_index=True) if all_data else pd.DataFrame()
    else:
        data = load_and_process_traces(file_pattern, condition_name)
    
    if not data.empty:
        condition_data[condition_name] = data
        condition_summaries[condition_name] = create_binned_average(data, bin_size=BIN_SIZE, position_range=PLOT_RANGE)
        print(f"  Loaded {data['sample'].nunique()} samples")
        print(f"  Generated {len(condition_summaries[condition_name])} bins")
    else:
        print(f"  No data found for {condition_name}")

# Create the publication-quality plot
fig, ax = plt.subplots(1, 1, figsize=(14, 7))

# Color palette matching thesis figure (p65 in gray, 2xTIMs-p65 in magenta/pink)
colors = ['#808080', '#e91e63']  # Gray for p65, magenta for 2xTIMs-p65

for i, (condition_name, summary) in enumerate(condition_summaries.items()):
    if summary.empty:
        continue
        
    color = colors[i % len(colors)]
    
    # Plot mean trace
    ax.plot(summary['bin'], summary['mean'], 
            linewidth=2.5, color=color, label=condition_name)
    
    # Add confidence interval band
    ax.fill_between(summary['bin'], summary['ci_lower'], summary['ci_upper'],
                    alpha=0.25, color=color)

# Add TSS and TES lines
ax.axvline(TSS, color='black', linestyle='--', alpha=0.8, linewidth=2)
ax.axvline(TES, color='black', linestyle='--', alpha=0.8, linewidth=2)

# Add TSS and TES labels
ax.text(TSS, ax.get_ylim()[1] * 0.95, 'TSS', ha='center', va='top', 
        fontweight='bold', fontsize=12, 
        bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
ax.text(TES, ax.get_ylim()[1] * 0.95, 'TES', ha='center', va='top', 
        fontweight='bold', fontsize=12,
        bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))

# Styling to match thesis figure
ax.set_xlabel('Position', fontsize=12, fontweight='bold')
ax.set_ylabel('Read Count', fontsize=12, fontweight='bold')
ax.set_title('PRO-seq Signal: p65 vs 2xTIMs-p65', fontsize=14, fontweight='bold')

# Set x-axis to show relative positions
x_ticks = np.arange(TSS - 500, TES + 501, 250)
x_labels = []
for x in x_ticks:
    if x == TSS:
        x_labels.append('TSS')
    elif x == TES:
        x_labels.append('TES')
    elif x < TSS:
        x_labels.append(f'-{abs(x-TSS)} bp')
    else:
        x_labels.append(f'+{x-TES} bp')

ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels, rotation=45)

# Styling
ax.set_xlim(PLOT_RANGE)
ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('p65_vs_2xTIMs_p65_proseq_traces.png', dpi=300, bbox_inches='tight')
plt.show()

print("\nPlot saved as: p65_vs_2xTIMs_p65_proseq_traces.png")

# Print summary statistics for each condition
print("\n=== CONDITION COMPARISON ===")
for condition_name, summary in condition_summaries.items():
    if summary.empty:
        continue
        
    print(f"\n{condition_name}:")
    print(f"  Samples: {condition_data[condition_name]['sample'].nunique()}")
    print(f"  Max signal: {summary['mean'].max():.1f}")
    print(f"  Max position: {summary.loc[summary['mean'].idxmax(), 'bin']:.0f}")
    
    # TSS region analysis
    tss_region = summary[(summary['bin'] >= TSS-50) & (summary['bin'] <= TSS+50)]
    if len(tss_region) > 0:
        tss_signal = tss_region['mean'].sum()
        print(f"  TSS signal (±50bp): {tss_signal:.1f}")

# Save all condition summaries
for condition_name, summary in condition_summaries.items():
    if not summary.empty:
        filename = f'{condition_name.lower().replace(" ", "_")}_summary.csv'
        summary.to_csv(filename, index=False)
        print(f"Saved {condition_name} summary to: {filename}")

print("\n=== INSTRUCTIONS FOR YOUR DATA ===")
print("1. Modify the 'conditions' dictionary to match your experimental design")
print("2. Group samples by treatment, timepoint, genotype, etc.")
print("3. Use meaningful condition names for the legend")
print("4. Adjust colors if you have more than 5 conditions")
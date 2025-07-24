#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
from scipy import stats
from scipy.ndimage import gaussian_filter1d
import seaborn as sns

# Set publication style
plt.style.use('default')
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'legend.fontsize': 12,
    'figure.figsize': (12, 6),
    'axes.linewidth': 1.5,
    'lines.linewidth': 2
})

def load_and_process_traces(file_pattern, condition_name):
    """Load trace files and return combined dataframe with normalization"""
    files = sorted(glob.glob(file_pattern))
    print(f"Found {len(files)} files for {condition_name}")
    
    all_data = []
    total_reads_per_sample = []
    
    for i, file in enumerate(files):
        df = pd.read_csv(file, sep='\t', header=None, names=['pos', 'count'])
        
        # Calculate total reads for this sample (for RPM normalization)
        total_reads = df['count'].sum()
        total_reads_per_sample.append(total_reads)
        
        # Normalize to Reads Per Million (RPM)
        df['count_rpm'] = (df['count'] / total_reads) * 1_000_000
        
        df['sample'] = f"{condition_name}_Rep{i+1}"
        df['condition'] = condition_name
        all_data.append(df)
    
    print(f"  Total reads per sample: {total_reads_per_sample}")
    return pd.concat(all_data, ignore_index=True) if all_data else pd.DataFrame()

def create_smooth_average(data, bin_size=10, position_range=None, smooth_sigma=2):
    """Create smoothed averages with confidence intervals"""
    if data.empty:
        return pd.DataFrame()
        
    if position_range:
        data = data[(data['pos'] >= position_range[0]) & (data['pos'] <= position_range[1])]
    
    # Create position grid
    min_pos = data['pos'].min()
    max_pos = data['pos'].max()
    
    # Align to bin_size
    bin_start = (min_pos // bin_size) * bin_size
    bin_end = ((max_pos // bin_size) + 1) * bin_size
    bins = range(int(bin_start), int(bin_end), bin_size)
    
    # Assign bins
    data['bin'] = pd.cut(data['pos'], bins=bins, labels=bins[:-1], include_lowest=True)
    data['bin'] = data['bin'].astype(float)
    
    # Calculate RPM statistics for each bin
    bin_stats = data.groupby(['bin', 'sample'])['count_rpm'].sum().reset_index()
    
    # Calculate mean and confidence intervals
    summary = bin_stats.groupby('bin')['count_rpm'].agg([
        'mean', 'std', 'count', 'sem'
    ]).reset_index()
    
    # Calculate 95% confidence intervals
    summary['ci_lower'] = summary['mean'] - 1.96 * summary['sem']
    summary['ci_upper'] = summary['mean'] + 1.96 * summary['sem']
    summary['ci_lower'] = np.maximum(summary['ci_lower'], 0)
    
    # Apply Gaussian smoothing to reduce noise (like in published figure)
    if smooth_sigma > 0:
        summary['mean_smooth'] = gaussian_filter1d(summary['mean'], sigma=smooth_sigma)
        summary['ci_lower_smooth'] = gaussian_filter1d(summary['ci_lower'], sigma=smooth_sigma)
        summary['ci_upper_smooth'] = gaussian_filter1d(summary['ci_upper'], sigma=smooth_sigma)
    else:
        summary['mean_smooth'] = summary['mean']
        summary['ci_lower_smooth'] = summary['ci_lower']
        summary['ci_upper_smooth'] = summary['ci_upper']
    
    return summary

# Parameters
TSS, TES = 8235, 8949
PLOT_RANGE = (TSS - 500, TES + 500)
BIN_SIZE = 10
SMOOTH_SIGMA = 2  # Gaussian smoothing parameter

print("=== Creating Publication-Style PRO-seq Plot ===")

# Load corrected data (strand-specific)
conditions = {
    'p65': 'traces/PR1567_S16_corrected.binned.tsv traces/PR1568_S17_corrected.binned.tsv',
    '2xTIMs-p65': 'traces/PR1569_S18_corrected.binned.tsv traces/PR1570_S19_corrected.binned.tsv'
}

# Process each condition
condition_summaries = {}
for condition_name, file_pattern in conditions.items():
    print(f"\nProcessing {condition_name}...")
    
    # Handle space-separated file lists
    files = file_pattern.split()
    all_data = []
    
    for i, file in enumerate(files):
        if glob.glob(file):
            df = pd.read_csv(file, sep='\t', header=None, names=['pos', 'count'])
            
            # Calculate total reads for RPM normalization
            total_reads = df['count'].sum()
            df['count_rpm'] = (df['count'] / total_reads) * 1_000_000
            
            df['sample'] = f"{condition_name}_Rep{i+1}"
            df['condition'] = condition_name
            all_data.append(df)
    
    if all_data:
        data = pd.concat(all_data, ignore_index=True)
        condition_summaries[condition_name] = create_smooth_average(
            data, bin_size=BIN_SIZE, position_range=PLOT_RANGE, smooth_sigma=SMOOTH_SIGMA
        )
        print(f"  Processed {len(all_data)} samples")

# Create publication-style plot
fig, ax = plt.subplots(1, 1, figsize=(14, 6))

# Colors matching the published figure
colors = {'p65': '#808080', '2xTIMs-p65': '#e91e63'}  # Gray and magenta

for condition_name, summary in condition_summaries.items():
    if summary.empty:
        continue
        
    color = colors[condition_name]
    
    # Plot smoothed mean trace
    ax.plot(summary['bin'], summary['mean_smooth'], 
            linewidth=2.5, color=color, label=condition_name, alpha=0.9)
    
    # Add smoothed confidence interval band
    ax.fill_between(summary['bin'], summary['ci_lower_smooth'], summary['ci_upper_smooth'],
                    alpha=0.25, color=color)

# Add TSS and TES lines
ax.axvline(TSS, color='black', linestyle='--', alpha=0.7, linewidth=1.5)
ax.axvline(TES, color='black', linestyle='--', alpha=0.7, linewidth=1.5)

# Add labels for TSS and TES
ax.text(TSS, ax.get_ylim()[1] * 0.95, 'TSS', ha='center', va='top', 
        fontweight='bold', fontsize=12)
ax.text(TES, ax.get_ylim()[1] * 0.95, 'TES', ha='center', va='top', 
        fontweight='bold', fontsize=12)

# Styling to match published figure
ax.set_xlabel('Position', fontsize=14, fontweight='bold')
ax.set_ylabel('Read Count (RPM)', fontsize=14, fontweight='bold')
ax.set_title('PRO-seq Signal: p65 vs 2xTIMs-p65', fontsize=16, fontweight='bold')

# Set x-axis to show relative positions like published figure
x_ticks = [TSS - 500, TSS, TES, TES + 500]
x_labels = ['-500 bp', 'TSS', 'TES', '+500 bp']
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels)

# Clean styling
ax.set_xlim(PLOT_RANGE)
ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=False, 
          title='Condition', title_fontsize=12)
ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set y-axis to start from 0
ax.set_ylim(bottom=0)

plt.tight_layout()
plt.savefig('publication_style_proseq_traces.png', dpi=300, bbox_inches='tight')
plt.show()

print("\nPublication-style plot saved as: publication_style_proseq_traces.png")

# Print normalization info
print("\n=== NORMALIZATION DETAILS ===")
for condition_name, file_pattern in conditions.items():
    files = file_pattern.split()
    print(f"\n{condition_name}:")
    for file in files:
        if glob.glob(file):
            df = pd.read_csv(file, sep='\t', header=None, names=['pos', 'count'])
            total_reads = df['count'].sum()
            print(f"  {file.split('/')[-1]}: {total_reads:,} total reads")

print("\n=== KEY IMPROVEMENTS ===")
print("1. ✅ RPM normalization (Reads Per Million)")
print("2. ✅ Gaussian smoothing (reduces noise)")
print("3. ✅ Confidence intervals across replicates")
print("4. ✅ Publication-style colors and formatting")
print("5. ✅ Cleaner axis labels and styling")
print("6. ✅ Uses corrected strand-specific data")
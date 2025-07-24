#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import glob
from scipy import stats

# Set publication style
plt.style.use('default')
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9
})

def load_sample_data(file_path, condition, replicate):
    """Load a single sample's trace data"""
    df = pd.read_csv(file_path, sep='\t', header=None, names=['pos', 'count'])
    total_reads = df['count'].sum()
    df['count_rpm'] = (df['count'] / total_reads) * 1_000_000
    df['condition'] = condition
    df['replicate'] = replicate
    return df

def calculate_regional_metrics(df, tss, tes):
    """Calculate average reads in different genomic regions"""
    regions = {
        'Prom. prox.': (tss - 50, tss + 150),
        'Gene body': (tss + 150, tes - 150),
        'Distal': (tes - 150, tes),
        'Post-term.': (tes, tes + 500)
    }
    
    regional_data = {}
    for region_name, (start, end) in regions.items():
        region_reads = df[(df['pos'] >= start) & (df['pos'] <= end)]
        
        if len(region_reads) > 0:
            avg_reads_rpm = region_reads['count_rpm'].mean()
        else:
            avg_reads_rpm = 0
        
        regional_data[region_name] = {'avg_reads_rpm': avg_reads_rpm}
    
    return regional_data

# Parameters
TSS, TES = 8235, 8949

# Load and process data
samples = [
    ('traces/PR1567_S16_corrected.binned.tsv', 'p65', 'Rep1'),
    ('traces/PR1568_S17_corrected.binned.tsv', 'p65', 'Rep2'),
    ('traces/PR1569_S18_corrected.binned.tsv', '2xTIMs-p65', 'Rep1'),
    ('traces/PR1570_S19_corrected.binned.tsv', '2xTIMs-p65', 'Rep2')
]

all_regional_data = []
for file_path, condition, replicate in samples:
    if glob.glob(file_path):
        df = load_sample_data(file_path, condition, replicate)
        regional_metrics = calculate_regional_metrics(df, TSS, TES)
        
        for region, metrics in regional_metrics.items():
            all_regional_data.append({
                'condition': condition,
                'replicate': replicate,
                'region': region,
                'avg_reads_rpm': metrics['avg_reads_rpm']
            })

results_df = pd.DataFrame(all_regional_data)

# Calculate normalized data
normalized_data = []
for condition in ['p65', '2xTIMs-p65']:
    for replicate in ['Rep1', 'Rep2']:
        sample_data = results_df[(results_df['condition'] == condition) & 
                                (results_df['replicate'] == replicate)]
        
        gene_body_signal = sample_data[sample_data['region'] == 'Gene body']['avg_reads_rpm'].iloc[0]
        
        if gene_body_signal > 0:
            for _, row in sample_data.iterrows():
                normalized_signal = row['avg_reads_rpm'] / gene_body_signal
                normalized_data.append({
                    'condition': condition,
                    'replicate': replicate,
                    'region': row['region'],
                    'normalized_signal': normalized_signal
                })

normalized_df = pd.DataFrame(normalized_data)

# Create complete figure
fig = plt.figure(figsize=(14, 10))

# Create gene diagram at top
ax_gene = plt.subplot2grid((3, 2), (0, 0), colspan=2, rowspan=1)

# Gene diagram parameters
gene_start = TSS - 500
gene_end = TES + 500
total_length = gene_end - gene_start

# Scale positions for plotting
def scale_pos(pos):
    return (pos - gene_start) / total_length

# Draw gene structure
gene_y = 0.5
gene_height = 0.3

# Main gene body (gray)
gene_rect = patches.Rectangle((scale_pos(TSS), gene_y - gene_height/2), 
                             scale_pos(TES) - scale_pos(TSS), gene_height,
                             facecolor='lightgray', edgecolor='black', linewidth=1)
ax_gene.add_patch(gene_rect)

# Add promoter region (light blue)
prom_rect = patches.Rectangle((scale_pos(TSS - 50), gene_y - gene_height/2), 
                             scale_pos(TSS + 150) - scale_pos(TSS - 50), gene_height,
                             facecolor='lightblue', edgecolor='black', linewidth=1, alpha=0.7)
ax_gene.add_patch(prom_rect)

# Add arrows for TSS and TES
arrow_props = dict(arrowstyle='->', lw=1.5, color='black')
ax_gene.annotate('', xy=(scale_pos(TSS), gene_y + gene_height/2 + 0.1), 
                xytext=(scale_pos(TSS), gene_y + gene_height/2 + 0.3),
                arrowprops=arrow_props)
ax_gene.annotate('', xy=(scale_pos(TES), gene_y + gene_height/2 + 0.1), 
                xytext=(scale_pos(TES), gene_y + gene_height/2 + 0.3),
                arrowprops=arrow_props)

# Add labels
ax_gene.text(scale_pos(TSS), gene_y + gene_height/2 + 0.4, 'TSS', ha='center', va='bottom', fontweight='bold')
ax_gene.text(scale_pos(TES), gene_y + gene_height/2 + 0.4, 'TES', ha='center', va='bottom', fontweight='bold')

# Add position markers
positions = [gene_start, TSS, TES, gene_end]
labels = ['-500 bp', 'TSS', 'TES', '+500 bp']
for pos, label in zip(positions, labels):
    ax_gene.axvline(scale_pos(pos), color='black', linestyle='--', alpha=0.5)
    ax_gene.text(scale_pos(pos), 0, label, ha='center', va='top', fontsize=9)

# Add region labels
regions_pos = {
    'Promoter proxy': (scale_pos(TSS - 50), scale_pos(TSS + 150)),
    'Gene Body': (scale_pos(TSS + 150), scale_pos(TES - 150)),
    'Distal': (scale_pos(TES - 150), scale_pos(TES)),
    'Post Term': (scale_pos(TES), scale_pos(TES + 500))
}

for region, (start, end) in regions_pos.items():
    mid_point = (start + end) / 2
    ax_gene.text(mid_point, gene_y - gene_height/2 - 0.2, region, 
                ha='center', va='top', fontsize=9, style='italic')

ax_gene.set_xlim(0, 1)
ax_gene.set_ylim(-0.5, 1)
ax_gene.set_title('Gene Structure and Analysis Regions', fontweight='bold', pad=20)
ax_gene.axis('off')

# Create data plots
colors = {'p65': '#808080', '2xTIMs-p65': '#e91e63'}
region_order = ['Prom. prox.', 'Gene body', 'Distal', 'Post-term.']

# Left panel: Raw signals
ax1 = plt.subplot2grid((3, 2), (1, 0), rowspan=2)

plot_data_raw = []
for condition in ['p65', '2xTIMs-p65']:
    for region in region_order:
        region_data = results_df[(results_df['condition'] == condition) & 
                                (results_df['region'] == region)]
        for _, row in region_data.iterrows():
            plot_data_raw.append({
                'condition': condition,
                'region': region,
                'value': row['avg_reads_rpm']
            })

raw_plot_df = pd.DataFrame(plot_data_raw)

# Plot with individual points and means
for i, region in enumerate(region_order):
    for j, condition in enumerate(['p65', '2xTIMs-p65']):
        data = raw_plot_df[(raw_plot_df['region'] == region) & 
                          (raw_plot_df['condition'] == condition)]
        
        x_pos = i + (j - 0.5) * 0.15
        
        # Individual points
        ax1.scatter([x_pos] * len(data), data['value'], 
                   color=colors[condition], alpha=0.8, s=80, 
                   label=condition if i == 0 else "", zorder=3)
        
        # Mean line
        mean_val = data['value'].mean()
        ax1.plot([x_pos - 0.08, x_pos + 0.08], [mean_val, mean_val], 
                color=colors[condition], linewidth=4, zorder=4)

ax1.set_xlabel('Genomic Region', fontweight='bold')
ax1.set_ylabel('avg window reads', fontweight='bold')
ax1.set_xticks(range(len(region_order)))
ax1.set_xticklabels(region_order, rotation=45, ha='right')
ax1.legend(title='Condition', loc='upper left')
ax1.grid(True, alpha=0.3, zorder=1)
ax1.set_ylim(bottom=0)

# Right panel: Normalized signals
ax2 = plt.subplot2grid((3, 2), (1, 1), rowspan=2)

plot_data_norm = []
for condition in ['p65', '2xTIMs-p65']:
    for region in region_order:
        region_data = normalized_df[(normalized_df['condition'] == condition) & 
                                   (normalized_df['region'] == region)]
        for _, row in region_data.iterrows():
            plot_data_norm.append({
                'condition': condition,
                'region': region,
                'value': row['normalized_signal']
            })

norm_plot_df = pd.DataFrame(plot_data_norm)

for i, region in enumerate(region_order):
    for j, condition in enumerate(['p65', '2xTIMs-p65']):
        data = norm_plot_df[(norm_plot_df['region'] == region) & 
                           (norm_plot_df['condition'] == condition)]
        
        x_pos = i + (j - 0.5) * 0.15
        
        # Individual points
        ax2.scatter([x_pos] * len(data), data['value'], 
                   color=colors[condition], alpha=0.8, s=80, zorder=3)
        
        # Mean line
        mean_val = data['value'].mean()
        ax2.plot([x_pos - 0.08, x_pos + 0.08], [mean_val, mean_val], 
                color=colors[condition], linewidth=4, zorder=4)

ax2.set_xlabel('Genomic Region', fontweight='bold')
ax2.set_ylabel('normalized avg window reads', fontweight='bold')
ax2.set_xticks(range(len(region_order)))
ax2.set_xticklabels(region_order, rotation=45, ha='right')
ax2.grid(True, alpha=0.3, zorder=1)
ax2.set_ylim(bottom=0)

plt.tight_layout()
plt.savefig('complete_regional_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

print("Complete regional analysis figure saved as: complete_regional_analysis.png")

# Print summary
print("\n=== SUMMARY STATISTICS ===")
for region in region_order:
    p65_raw = raw_plot_df[(raw_plot_df['region'] == region) & 
                         (raw_plot_df['condition'] == 'p65')]['value']
    tims_raw = raw_plot_df[(raw_plot_df['region'] == region) & 
                          (raw_plot_df['condition'] == '2xTIMs-p65')]['value']
    
    if len(p65_raw) > 0 and len(tims_raw) > 0:
        fold_change = tims_raw.mean() / p65_raw.mean() if p65_raw.mean() > 0 else float('inf')
        print(f"{region}: {fold_change:.2f}-fold increase with 2xTIMs-p65")

results_df.to_csv('regional_analysis_detailed.csv', index=False)
print(f"\nResults saved to regional_analysis_detailed.csv")
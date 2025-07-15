import pyBigWig
import numpy as np
import matplotlib.pyplot as plt

# Input BigWigs — update paths if needed
bigwigs = {
    "p65": [
        "results/bigwig/PR1567_S16_R1.bw",
        "results/bigwig/PR1569_S18_R1.bw"
    ],
    "2xTIMs-p65": [
        "results/bigwig/PR1568_S17_R1.bw",
        "results/bigwig/PR1570_S19_R1.bw"
    ]
}

# Function to auto-detect plasmid contig
def detect_plasmid_contig(bw):
    chroms = bw.chroms()
    for chrom in chroms:
        if not chrom.startswith('chr') and not chrom.startswith('KI') and not chrom.startswith('GL'):
            return chrom, chroms[chrom]
    return None, None

# Detect plasmid contig
example_bw = pyBigWig.open(next(iter(bigwigs.values()))[0])
plasmid_chrom, plasmid_len = detect_plasmid_contig(example_bw)
example_bw.close()

if plasmid_chrom is None:
    raise RuntimeError("Plasmid contig not detected! Please check your BigWig files.")

print(f"Detected plasmid contig: {plasmid_chrom} (length: {plasmid_len})")

# Plot setup
plt.figure(figsize=(14, 6))
window = 50  # smoothing window
start_pos = 6000

# Colors for the lines
colors = {'p65': '#1f77b4', '2xTIMs-p65': '#ff7f0e'}

for label, paths in bigwigs.items():
    all_values = []
    
    # Load all replicates
    for bw_path in paths:
        bw = pyBigWig.open(bw_path)
        vals = np.array(bw.values(plasmid_chrom, 0, plasmid_len))
        bw.close()
        all_values.append(vals)
    
    all_values = np.vstack(all_values)
    
    # Smooth each replicate individually
    smoothed_reps = []
    for vals in all_values:
        smoothed = np.convolve(vals, np.ones(window)/window, mode='same')
        smoothed_reps.append(smoothed)
    
    smoothed_reps = np.vstack(smoothed_reps)
    
    # Extract region from 6000bp to end
    region_data = smoothed_reps[:, start_pos:]
    region_x = np.arange(start_pos, plasmid_len)
    
    # Calculate mean and standard deviation
    mean_vals = np.mean(region_data, axis=0)
    std_vals = np.std(region_data, axis=0)
    
    # Plot mean line
    color = colors[label]
    plt.plot(region_x, mean_vals, label=label, color=color, linewidth=2)
    
    # Add shaded standard deviation
    plt.fill_between(region_x, mean_vals - std_vals, mean_vals + std_vals, 
                     color=color, alpha=0.2)

# TSS and TES positions
TSS = 8236
TES = 8949

# Add vertical lines for TSS and TES with clear labels
plt.axvline(x=TSS, color='gray', linestyle='--', linewidth=1.5, alpha=0.8)
plt.axvline(x=TES, color='black', linestyle='--', linewidth=1.5, alpha=0.8)

# Add text labels for TSS and TES
plt.text(TSS, plt.ylim()[1] * 0.95, 'mCherry TSS', rotation=90, 
         verticalalignment='top', horizontalalignment='right', 
         fontsize=10, color='gray', fontweight='bold')
plt.text(TES, plt.ylim()[1] * 0.95, 'TES/Puro start', rotation=90, 
         verticalalignment='top', horizontalalignment='right', 
         fontsize=10, color='black', fontweight='bold')

# Add shaded region for mCherry
plt.axvspan(TSS, TES, color='lightgray', alpha=0.3, label='mCherry')

# Axis labels and title
plt.xlabel("Plasmid Position (bp)", fontsize=12)
plt.ylabel("Coverage (CPM)", fontsize=12)
plt.title("PRO-seq Signal Across Plasmid Cassette", fontsize=14, fontweight='bold')

# Set x-axis limits to show only the region of interest
plt.xlim(start_pos, plasmid_len)

# Customize x-axis ticks
x_ticks = np.arange(6000, plasmid_len + 1000, 1000)
x_ticks = x_ticks[x_ticks <= plasmid_len]
plt.xticks(x_ticks)

# Add grid for better readability
plt.grid(True, alpha=0.3)

# Legend
plt.legend(loc="upper left", frameon=True, fancybox=True, shadow=True)

# Layout and save
plt.tight_layout()
plt.savefig("plasmid_cassette_signal.png", dpi=300, bbox_inches='tight')


print("✅ Plot saved as plasmid_cassette_signal.png")
print(f"✅ Plot shows region from {start_pos}bp to {plasmid_len}bp with standard deviation shading")
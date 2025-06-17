import pyBigWig
import numpy as np
import matplotlib.pyplot as plt

# Input BigWigs
bigwigs = {
    "p65": "results/bigwig/PR1567_S16_R1.bw",
    "2xTIMs-p65": "results/bigwig/PR1568_S17_R1.bw"
}

# Function to auto-detect plasmid contig (assuming non-chromosomal name)
def detect_plasmid_contig(bw):
    chroms = bw.chroms()
    # Look for contig that doesn't match standard chromosomes
    for chrom in chroms:
        if not chrom.startswith('chr') and not chrom.startswith('KI') and not chrom.startswith('GL'):
            return chrom, chroms[chrom]
    return None, None

# Open first BigWig to detect plasmid contig
example_bw = pyBigWig.open(next(iter(bigwigs.values())))
plasmid_chrom, plasmid_len = detect_plasmid_contig(example_bw)
example_bw.close()

if plasmid_chrom is None:
    raise RuntimeError("Plasmid contig not detected! Please check your BigWig files.")

print(f"Detected plasmid contig: {plasmid_chrom} (length: {plasmid_len})")

# Prepare plotting
x = np.arange(0, plasmid_len)
plt.figure(figsize=(10,4))

for label, bw_path in bigwigs.items():
    bw = pyBigWig.open(bw_path)
    values = np.array(bw.values(plasmid_chrom, 0, plasmid_len))
    bw.close()

    # Optional smoothing
    window = 50
    smoothed = np.convolve(values, np.ones(window)/window, mode='same')

    plt.plot(x, smoothed, label=label)

plt.xlabel("Plasmid Position (bp)")
plt.ylabel("Coverage (CPM)")
plt.title("PRO-seq Signal Across Plasmid")
plt.legend()
plt.tight_layout()

# Save to file
plt.savefig("plasmid_proseq_signal.png", dpi=300)
print("✅ Plot saved as plasmid_proseq_signal.png")


import pyBigWig
import numpy as np
import matplotlib.pyplot as plt

# Setup
bigwigs = {
    "p65": "results/bigwig/PR1567_S16_R1.bw",
    "2xTIMs-p65": "results/bigwig/PR1568_S17_R1.bw"
}

chrom = "ProSEQ_project/refs/psc146-8xzf10bs-ybtata-mche-hpgk-puro.fasta"  # Adjust to match your reference
region_start = 0
region_end = 10707  # or whatever covers your plasmid

x = np.arange(region_start, region_end)

plt.figure(figsize=(10,4))

for label, bw_path in bigwigs.items():
    bw = pyBigWig.open(bw_path)
    values = np.array(bw.values(chrom, region_start, region_end))
    bw.close()

    # Optional: apply smoothing
    window = 50
    smoothed = np.convolve(values, np.ones(window)/window, mode='same')

    plt.plot(x, smoothed, label=label)

plt.xlabel("Plasmid Position (bp)")
plt.ylabel("Coverage (CPM)")
plt.title("PRO-seq Signal Across Plasmid")
plt.legend()
plt.tight_layout()
plt.show()

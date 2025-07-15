#!/usr/bin/env python3
"""
Plot PRO-seq R2 signal across a plasmid cassette and display a feature map.

• Uses the *R2* bigWig files listed in `bigwigs`.
• Automatically detects the plasmid contig in each bigWig (assumes it is the
  only chromosome that does **not** start with “chr”, “KI”, or “GL”).
• Smooths each replicate with a sliding-window mean (`window` bp).
• Plots the mean ± SD of replicates for each condition.
• Adds an inset axis underneath showing major plasmid features (TSS, TES,
  mCherry ORF, promoter, PuroR, …).

Requires:  pyBigWig, numpy, matplotlib  (pip/conda install as needed).
"""

# ------------------------------------------------------------------ imports ----
import os
import sys
import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ------------------------------------------------------------------- inputs ----
bigwigs = {
    "p65": [
        "results/bigwig/PR1567_S16_R2.bw",
        "results/bigwig/PR1569_S18_R2.bw",
    ],
    "2xTIMs-p65": [
        "results/bigwig/PR1568_S17_R2.bw",
        "results/bigwig/PR1570_S19_R2.bw",
    ],
}

# parameters you might want to tweak
window     = 50         # smoothing window (bp)
start_pos  = 6000       # left edge of region to plot
colors     = {"p65": "#1f77b4", "2xTIMs-p65": "#ff7f0e"}

# major genomic landmarks / features on the plasmid
features = [
    (8236, 8236,  "TSS",        "grey"),       # single-bp feature → vertical line
    (8949, 8949,  "TES",        "black"),
    (8239, 8943,  "mCherry",    "#b90746"),
    (8969, 9469,  "hPGK prom.", "#c7b0e3"),
    (9479, 10078, "PuroR",      "#b4abac"),
    (7774, 8158,  "8x ZF10", "red"),
    (10116, 10707, "WPRE", "blue"),
    
]

# ------------------------------------------------------------------ helpers ----
def detect_plasmid_contig(bw: pyBigWig.pyBigWig):
    """Return (contig_name, length) of the plasmid in a bigWig file."""
    for chrom, length in bw.chroms().items():
        if not chrom.startswith(("chr", "KI", "GL")):
            return chrom, length
    return None, None


def load_and_smooth(bw_path: str, chrom: str, length: int, win: int):
    """Load coverage vector from `bw_path`, smooth with `win` bp mean filter."""
    with pyBigWig.open(bw_path) as bw:
        vals = np.array(bw.values(chrom, 0, length), dtype=float)
    vals = np.nan_to_num(vals, nan=0.0)               # replace NaNs with 0
    return np.convolve(vals, np.ones(win) / win, mode="same")


# ---------------------------------------------------------------- detect contig
# Use the first bigWig we find to auto-detect the plasmid chromosome name
try:
    first_bw = next(iter(next(iter(bigwigs.values()))))
except StopIteration:
    sys.exit("❌ No bigWig paths supplied!")

with pyBigWig.open(first_bw) as bw:
    plasmid_chrom, plasmid_len = detect_plasmid_contig(bw)

if plasmid_chrom is None:
    sys.exit("❌ Plasmid contig not detected in the bigWig(s).")
print(f"Detected plasmid contig: {plasmid_chrom}  (length = {plasmid_len:,} bp)")

# ------------------------------------------------------------ figure & axes ----
fig, ax_cov = plt.subplots(figsize=(14, 6))  # main coverage axis

# ------------------------------------------------------- plot coverage data ----
for label, bw_paths in bigwigs.items():
    smoothed_reps = []

    for path in bw_paths:
        if not os.path.exists(path):
            print(f"⚠️  Skipping missing file: {path}", file=sys.stderr)
            continue
        smoothed = load_and_smooth(path, plasmid_chrom, plasmid_len, window)
        smoothed_reps.append(smoothed)

    if not smoothed_reps:
        print(f"⚠️  No data for condition '{label}'.", file=sys.stderr)
        continue

    smoothed_reps = np.vstack(smoothed_reps)

    # restrict to region of interest
    region_x     = np.arange(start_pos, plasmid_len)
    region_data  = smoothed_reps[:, start_pos:]
    mean_vals    = region_data.mean(axis=0)
    std_vals     = region_data.std(axis=0)

    col = colors.get(label, None)
    ax_cov.plot(region_x, mean_vals, label=label, color=col, lw=2)
    ax_cov.fill_between(region_x, mean_vals - std_vals, mean_vals + std_vals,
                        color=col, alpha=0.25)

# --------------------------------------------------- add feature track axis ----
ax_feat = ax_cov.inset_axes([0, -0.14, 1, 0.08], sharex=ax_cov)
ax_feat.set_ylim(0, 1)
ax_feat.axis("off")

for s, e, lbl, col in features:
    if s == e:  # single-bp feature (TSS/TES)
        ax_feat.axvline(s, color=col, lw=2)
        ax_feat.text(s, 0.6, lbl, ha="center", va="bottom", fontsize=8, rotation=90)
    else:       # block feature (genes, promoters, etc.)
        ax_feat.add_patch(
            mpatches.Rectangle((s, 0.2), e - s, 0.6, facecolor=col, edgecolor="none", alpha=0.6)
        )
        ax_feat.text(s + (e - s) / 2, 0.1, lbl, ha="center", va="top", fontsize=8)

# --------------------------------------------- annotate TSS / TES on coverage ----
TSS, TES = 8236, 8949
ax_cov.axvline(TSS, color="gray",  ls="--", lw=1.5, alpha=0.8)
ax_cov.axvline(TES, color="black", ls="--", lw=1.5, alpha=0.8)
ax_cov.text(TSS, ax_cov.get_ylim()[1] * 0.95, "mCherry TSS",
            rotation=90, ha="right", va="top", fontsize=10, color="gray",  fontweight="bold")
ax_cov.text(TES, ax_cov.get_ylim()[1] * 0.95, "TES / Puro start",
            rotation=90, ha="right", va="top", fontsize=10, color="black", fontweight="bold")
ax_cov.axvspan(TSS, TES, color="lightgray", alpha=0.25, label="mCherry")

# -------------------------------------------------------- prettify & export ----
ax_cov.set_xlabel("Plasmid Position (bp)", fontsize=12)
ax_cov.set_ylabel("Coverage (CPM)",        fontsize=12)
ax_cov.set_title("PRO-seq R2 Signal Across Plasmid Cassette",
                 fontsize=14, fontweight="bold")
ax_cov.set_xlim(start_pos, plasmid_len)
x_ticks = np.arange(6000, plasmid_len + 1000, 1000)
ax_cov.set_xticks(x_ticks)
ax_cov.grid(True, alpha=0.3)
ax_cov.legend(loc="upper left", frameon=True, fancybox=True, shadow=True)

plt.tight_layout()
fig.subplots_adjust(bottom=0.18)     # room for the feature map
out_png = "plasmid_cassette_R2_signal.png"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"✅ Plot with feature map saved ➜  {out_png}")

#!/usr/bin/env python3
"""
Produce 10-nt binned 3′ trace, PI, and window means for a PRO-seq library.

usage: proseq_count.py <in.bed> <TSS> <TES> <BIN> <prefix>
"""
import sys, pandas as pd, numpy as np, pathlib

bed, TSS, TES, BIN, prefix = sys.argv[1:]
TSS, TES, BIN = map(int, (TSS, TES, BIN))

df = pd.read_csv(bed, sep='\t', header=None,
                 names=['chr','start','end','name','score','strand'])

# Infer 3′ coordinate for each alignment ---------------------
df['threep'] = np.where(df.strand == '-', df.end - 1, df.start)

# Raw positional counts --------------------------------------
counts = df.threep.value_counts().sort_index()

# 10-nt binning ----------------------------------------------
binned = (counts.groupby((counts.index // BIN) * BIN).sum()
                 .rename('count').reset_index()
                 .rename(columns={'index':'bin_start'}))
binned.to_csv(f"{prefix}.binned_counts.tsv", sep='\t', index=False)

# Pausing index ----------------------------------------------
pi_numer = counts.loc[(TSS-50):(TSS+50)].sum()
pi_denom = counts.loc[(TSS+50):(TSS+250)].sum()
pi = pi_numer / pi_denom if pi_denom else np.nan

# Window means and normalisation -----------------------------
windows = {
    "promoter"  : (TSS,       TSS+150),
    "gene_body" : (TSS+150,   TES-150),
    "pre_TES"   : (TES-150,   TES),
    "post_TES"  : (TES,       TES+150),
    "downstream": (TES+150,   TES+500)
}

means = {k: counts.loc[v[0]:v[1]].sum() / (v[1]-v[0]) for k,v in windows.items()}
norm_factor = means["gene_body"] or np.nan
norm_means = {f"{k}_norm": v / norm_factor for k,v in means.items()}

with open(f"{prefix}.metrics.txt",'w') as fh:
    fh.write(f"PI\t{pi}\n")
    for k,v in means.items():
        fh.write(f"{k}_mean\t{v}\n")
    for k,v in norm_means.items():
        fh.write(f"{k}\t{v}\n")

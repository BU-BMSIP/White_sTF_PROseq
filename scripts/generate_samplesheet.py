#!/usr/bin/env python3

import os
import csv
from pathlib import Path

# Set paths
data_dir = Path("data")
output_csv = Path("full_samplesheet.csv")

# Gather all FASTQ files
fastq_files = sorted(data_dir.glob("*.fastq.gz"))

# Parse into rows
rows = []
for fq in fastq_files:
    # Extract sample name from filename
    # Example: LIB063021_SRN00286021_Khalil001-PR1567-656onBR1_S16_R1_001.fastq.gz
    # Suggested sample name: PR1567-656onBR1_S16_R1
    parts = fq.name.split("_")
    if len(parts) < 5:
        print(f"⚠️ Skipping unrecognized filename: {fq.name}")
        continue
    sample_name = f"{parts[2].split('-')[1]}_{parts[3]}_{parts[4]}"  # e.g. PR1567_S16_R1
    rows.append({"name": sample_name, "path": str(fq)})

# Write CSV
with open(output_csv, "w", newline="") as csvfile:
    fieldnames = ["name", "path"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"✅ Wrote {len(rows)} entries to {output_csv}")

#!/usr/bin/env bash

echo "=== Finding Sample 19 ==="

# First, let's search for any files related to sample 19
echo "Searching for PR1570_S19 files..."
find . -name "*PR1570_S19*" -type f 2>/dev/null | head -20

echo ""
echo "Searching for any S19 files..."
find . -name "*S19*" -type f 2>/dev/null | head -20

echo ""
echo "Checking data directory contents..."
ls -la data/ | grep -E "(S19|19)"

echo ""
echo "Checking if there are any numbered samples around 19..."
ls -la data/ | grep -E "S1[89]|S2[01]" | head -10

echo ""
echo "=== Current BED files available ==="
ls -la bed/ | head -20

echo ""
echo "=== Checking for any alignment files that might contain S19 ==="
find . -name "*S19*" -o -name "*19*" | grep -v ".nextflow" | head -20

# Let's also check what samples we actually have
echo ""
echo "=== Available samples in bed/ directory ==="
ls bed/*.bed 2>/dev/null | sed 's/.*\///; s/\.bed$//' | sort

echo ""
echo "=== Available samples in data/ directory ==="
ls data/*.fastq.gz 2>/dev/null | sed 's/.*\///; s/\.fastq\.gz$//' | sort

echo ""
echo "=== Let's see if there's a pattern in the sample numbering ==="
ls data/ | grep -E "S[0-9]+" | sed 's/.*S\([0-9]\+\).*/\1/' | sort -n | tail -10

echo ""
echo "=== If sample 19 exists in data/ but not in bed/, we need to process it ==="
if [ -f "data/PR1570_S19.fastq.gz" ]; then
    echo "Found PR1570_S19.fastq.gz in data/ directory!"
    echo "File size: $(ls -lh data/PR1570_S19.fastq.gz | awk '{print $5}')"
    echo "This sample needs to be processed through the alignment pipeline."
    echo ""
    echo "To process this sample, you would need to:"
    echo "1. Run the alignment pipeline on PR1570_S19.fastq.gz"
    echo "2. Generate the corresponding BED file"
    echo "3. Then run the trace generation script"
elif [ -f "data/PR1570_S19_R1.fastq.gz" ] || [ -f "data/PR1570_S19_001.fastq.gz" ]; then
    echo "Found PR1570_S19 with different naming convention!"
    ls -la data/*S19* 2>/dev/null
else
    echo "Sample 19 not found in data/ directory"
    echo "Checking if it might be named differently..."
    ls -la data/ | grep -i "19"
fi
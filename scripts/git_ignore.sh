# Create comprehensive .gitignore
cat > .gitignore << 'EOF'
# Nextflow cache and logs
.nextflow/
.nextflow.log
work/

# Large data files
data/
trimmed/
bam/
dedup/
umi/
results/

# File types to ignore
*.bam
*.bai
*.fastq
*.fastq.gz
*.fq.gz
*.sam

# Temporary files
unmapped.*
q
q.pub

# OS files
.DS_Store
Thumbs.db

# Editor files
*.swp
*.swo
*~

# Spike-in intermediate files (these are large)
spike_in_counts/
EOF
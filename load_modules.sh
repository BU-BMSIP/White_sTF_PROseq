# load_everything.sh
module purge
module load miniconda/23.11.0
module load snakemake/8.11.3
mamba activate /share/pkg.8/snakemake/8.11.3/install/snakemake
module load salmon
module load fastqc
module load star
module load subread

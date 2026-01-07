#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J Swe_DESeq2
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL

module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
source ../../.env

dir="$DIR/results/06_de/01_swe"
cd $dir
matrix_file="salmon.gene.counts.matrix"
col_data_file="coldata.tsv"

Rscript $SLURM_SUBMIT_DIR/diff_expression_swe.R $matrix_file $col_data_file
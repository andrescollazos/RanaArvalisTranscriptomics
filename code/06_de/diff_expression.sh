#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J DESeq2_Array
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH --array=0-5
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL

module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
source ../../.env

matrix_file="$DIR/data/quantification/salmon.gene.counts.matrix"
col_data_file="$DIR/data/diff_expression/coldata.tsv"
analyses=(
    "01_Temp_Swe"
    "02_Temp_Non_Swe"
    "03_NvS_Swe"
    "04_NvS_All"
    "05_Region"
    "06_Temp_Region"
)
analysis="${analyses[$SLURM_ARRAY_TASK_ID]}"
dir="$DIR/results/06_de/$analysis"
cd $dir

Rscript "${SLURM_SUBMIT_DIR}/${analysis}.R" $matrix_file $col_data_file
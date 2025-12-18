#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J PCAs
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 0:10:00
#SBATCH -p node
#SBATCH --array=0-15
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools R/4.3.1 R_packages/4.3.1
source ../../.env

cd "$DIR/results/05_quality_check/multi_qc_251215"
CONFIG="04_multi_qc.config.tsv"
read analysis < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$CONFIG" | cut -f1)

samples_config_file="samples.config.txt"
matrix_file="$analysis/${analysis}.centered.dat"

echo "Running Rscript for $analysis"
Rscript $SLURM_SUBMIT_DIR/05_pca_script.R $matrix_file $samples_config_file $analysis
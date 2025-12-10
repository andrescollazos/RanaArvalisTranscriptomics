#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J PtRHistogram
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 4:00:00
#SBATCH -p node
#SBATCH --array=0-1
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

# Define arrays
matrices=(
    "$DIR/data/quantification/salmon.gene.counts.matrix"
    "$DIR/data/quantification/salmon.isoform.counts.matrix"
)

folders=(
    "gene"
    "transcript"
)

out="$DIR/results/05_quality_check/all_samples"
samples="$out/full_samples.txt"
cd $out

# Pick the matrix and folder based on the array index
matrix=${matrices[$SLURM_ARRAY_TASK_ID]}
folder=${folders[$SLURM_ARRAY_TASK_ID]}

mkdir -p $folder
cd $folder

# Run the PtR command (same for both levels)
singularity exec --cleanenv \
    --env LC_ALL=C \
    $TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
    /usr/local/bin/Analysis/DifferentialExpression/PtR \
        --matrix $matrix \
        --samples $samples \
        --log2 --CPM --min_rowSums 10 \
        --barplot_sum_counts
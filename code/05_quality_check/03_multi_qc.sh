#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J PtRMultiQC
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 8:00:00
#SBATCH -p node
#SBATCH --array=0-13
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

matrix="$DIR/data/quantification/salmon.gene.counts.matrix"
dir="$DIR/results/05_quality_check/multi_qc"
analyses=(
    "swedish-only"
	"swedish-popxtemp"
    "non-swedish"
    "north-south"
    "temp-15-all_pop"
    "temp-20-all_pop"
    "pop-E-temp"
    "pop-C_Fin-temp"
    "pop-L-temp"
    "pop-Ka-temp"
    "pop-NA-temp"
    "pop-NL-temp"
    "pop-Upp-temp"
    "pop-VF-temp"
)
analysis="${analyses[$SLURM_ARRAY_TASK_ID]}"
out="$dir/$analysis"
cd $out

samples="${analysis}.samples.txt"

echo "Compare replicates"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $matrix \
		--samples $samples \
		--log2 --CPM \
		--min_rowSums 10 \
		--compare_replicates \
		--output $analysis

echo "Sample correlation matrix"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $matrix \
		--samples $samples \
		--log2 --CPM \
		--min_rowSums 10 \
		--sample_cor_matrix \
		--output $analysis

echo "Principal Component Analysis"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $matrix \
		--samples $samples \
		--log2 --CPM \
		--min_rowSums 10 \
		--center_rows --prin_comp 3 \
		--output $analysis
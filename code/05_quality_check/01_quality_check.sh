#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J PtRQualityCheck
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 4:00:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

matrix="$DIR/data/quantification/salmon.gene.counts.matrix"
out="$DIR/results/05_quality_check"
cd $out

singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $matrix \
		--samples full_samples.txt \
		--min_rowSums 10 \
		--compare_replicates
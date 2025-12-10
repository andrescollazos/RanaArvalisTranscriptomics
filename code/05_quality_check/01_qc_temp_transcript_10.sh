#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J 10_Temp_Transcript_PtRQC
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 8:00:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

transcript_matrix="$DIR/data/quantification/salmon.isoform.counts.matrix"
out="$DIR/results/05_quality_check/temperature"
samples="$out/full_samples.temperature.txt"
cd $out

#-----------------------------------------------------------------------------------
echo "TRANSCRIPT LEVEL:"
mkdir -p 10_transcript
cd 10_transcript

echo "Compare replicates"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $transcript_matrix \
		--samples $samples \
		--log2 --CPM \
		--min_rowSums 10 \
		--min_gene_expr_val 10 \
		--compare_replicates

echo "Sample correlation matrix"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $transcript_matrix \
		--samples $samples \
		--log2 --CPM \
		--min_rowSums 10 \
		--min_gene_expr_val 10 \
		--sample_cor_matrix

echo "Principal Component Analysis"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $transcript_matrix \
		--samples $samples \
		--log2 --CPM \
		--min_rowSums 10 \
		--min_gene_expr_val 10 \
		--center_rows --prin_comp 3
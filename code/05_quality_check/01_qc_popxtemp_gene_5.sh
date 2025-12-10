#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J 5_PopxTemp_Gene_PtRQC
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 8:00:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

gene_matrix="$DIR/data/quantification/salmon.gene.counts.matrix"
out="$DIR/results/05_quality_check/populationxtemperature"
samples="$out/full_samples.poptemp.txt"
cd $out

#-----------------------------------------------------------------------------------
echo "GENE LEVEL:"
mkdir -p 05_gene
cd 05_gene

echo "Compare replicates"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $gene_matrix \
		--samples $samples \
		--log2 --CPM \
		--min_rowSums 10 \
		--min_gene_expr_val 5 \
		--compare_replicates

echo "Sample correlation matrix"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $gene_matrix \
		--samples $samples \
		--log2 --CPM \
		--min_rowSums 10 \
		--min_gene_expr_val 5 \
		--sample_cor_matrix

echo "Principal Component Analysis"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $gene_matrix \
		--samples $samples \
		--log2 --CPM \
		--min_rowSums 10 \
		--min_gene_expr_val 5 \
		--prin_comp 3
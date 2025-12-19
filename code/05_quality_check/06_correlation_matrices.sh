#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J PtRQCBaltic
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 4:00:00
#SBATCH -p node
#SBATCH --array=0-7
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

matrix="$DIR/data/quantification/salmon.gene.counts.matrix"
dir="$DIR/results/05_quality_check/multi_qc_251215"
cd $dir
CONFIG="04_multi_qc.config.tsv"
read analysis samples_type min_rowSums top_genes top_variable_genes var_gene_method < <(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" $CONFIG)

cp baltic.samples.txt $analysis
cd $analysis

echo "Sample correlation matrix"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/Analysis/DifferentialExpression/PtR --matrix $matrix \
		--samples baltic.samples.txt \
		--log2 --CPM \
		--min_rowSums $min_rowSums \
		--top_genes $top_genes \
		--top_variable_genes $top_variable_genes \
		--var_gene_method $var_gene_method \
		--sample_cor_matrix \
		--output "${analysis}.baltic.cor_matrix"
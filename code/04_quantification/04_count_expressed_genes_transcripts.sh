#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J CountExpGeneTranscripts
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 1:00:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

out="$DIR/results/04_quantification"
cd $out

echo "Counting Numbers of Expressed Genes"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
        salmon.gene.TPM.not_cross_norm | tee salmon.gene.TPM.not_cross_norm.counts_by_min_TPM

echo "Counting Numbers of Expressed Transcripts"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
        salmon.isoform.TPM.not_cross_norm | tee salmon.isoform.TPM.not_cross_norm.counts_by_min_TPM
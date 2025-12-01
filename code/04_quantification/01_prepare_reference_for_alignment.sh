#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J PrepareReference
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 0:10:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

out="$DIR/results/04_quantification"
cd $out

singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/util/align_and_estimate_abundance.pl \
		--transcripts Trinity.fasta \
		--est_method salmon \
		--SS_lib_type FR \
		--thread_count $SLURM_CPUS_PER_TASK \
		--gene_trans_map Trinity.fasta.gene_trans_map \
		--prep_reference
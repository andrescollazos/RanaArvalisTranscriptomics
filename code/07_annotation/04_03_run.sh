#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J TrinotateRun
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 33:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL

source ../../.env

singularity exec \
	--bind $DIR:$DIR:rw \
	--cleanenv \
	--env LC_ALL=C \
	$TRINOTATE_SINGULARITY/trinotate.v4.0.2.simg \
	/usr/local/src/Trinotate/Trinotate \
	--db RanaArvalis_Trinotate.sqlite \
	--CPU $SLURM_CPUS_PER_TASK \
	--transcript_fasta $DIR/results/07_annotation/Trinity.fasta \
	--transdecoder_pep $DIR/results/07_annotation/Trinity.fasta.TD2.pep \
	--trinotate_data_dir $DIR/results/07_annotation/trinotate_data \
	--run "swissprot_blastp swissprot_blastx pfam infernal" \
	--use_diamond
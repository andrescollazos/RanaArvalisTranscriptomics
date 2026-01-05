#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J TrinotateInitDB
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL

source ../../.env

singularity exec \
	--bind $DIR:$DIR:rw \
	--cleanenv \
	--env LC_ALL=C \
	$TRINOTATE_SINGULARITY/trinotate.v4.0.2.simg \
	/usr/local/src/Trinotate/Trinotate --db RanaArvalis_Trinotate.sqlite \
	--init \
	--gene_trans_map $DIR/results/07_annotation/Trinity.fasta.gene_trans_map \
    --transcript_fasta $DIR/results/07_annotation/Trinity.fasta \
	--transdecoder_pep $DIR/results/07_annotation/Trinity.fasta.TD2.pep
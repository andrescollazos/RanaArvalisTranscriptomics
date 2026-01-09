#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J TrinotateReport
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 4:00:00
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
	--report > RanaArvalis_Trinotate.annotation.tsv

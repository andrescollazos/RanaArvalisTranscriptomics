#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J TrinotateCreateDB
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL

source ../../.env

singularity exec \
	--bind Build_Trinotate_Boilerplate_SQLite_db.pl:/usr/local/src/Trinotate/util/admin/Build_Trinotate_Boilerplate_SQLite_db.pl \
	--bind $DIR:$DIR:rw \
	--cleanenv \
	--env LC_ALL=C \
	$TRINOTATE_SINGULARITY/trinotate.v4.0.2.simg \
	/usr/local/src/Trinotate/Trinotate --create \
	--db RanaArvalis_Trinotate.sqlite \
	--trinotate_data_dir $DIR/results/07_annotation/trinotate_data \
	--use_diamond
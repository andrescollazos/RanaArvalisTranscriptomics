#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Contig_Nx_Stats
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 1:30:00
#SBATCH -p core
#SBATCH --cpus-per-task=2
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

dir="$DIR/results/03_evaluation/01_contig_nx_stats"
cd $dir

singularity exec --cleanenv \
    --env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/util/TrinityStats.pl Trinity.fasta
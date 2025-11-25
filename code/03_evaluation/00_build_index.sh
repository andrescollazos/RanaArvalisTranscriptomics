#!/bin/bash
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Eval_Read_Rep_Transcriptome
#SBATCH -t 4:00:00
#SBATCH -p node
#SBATCH --cpus-per-task=16
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools bowtie2/2.5.2 samtools/1.20
source ../../.env

transcriptome="$DIR/data/transcriptome/Trinity.fasta"
out="$DIR/results/03_evaluation/00_read_representation"

# Build a bowtie2 index for the transcriptome:
cd $out
bowtie2-build $transcriptome Trinity.fasta --threads $SLURM_CPUS_PER_TASK

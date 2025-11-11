#!/bin/bash
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J FastQC
#SBATCH -e %x.er
#SBATCH -t 5:00:00
#SBATCH -n 8
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load FastQC/0.11.9

source ../../.env
OUT="$DIR/results/01_preprocessing/fastqc_raw"
READS="$DIR/data/raw_reads/P32262_*/*.fastq.gz"

fastqc $READS -o $OUT -t 8

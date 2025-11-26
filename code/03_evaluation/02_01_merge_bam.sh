#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Merge_Sorted_Test
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 3:00:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools samtools/1.20
source ../../.env

out="$DIR/results/03_evaluation/00_read_representation"
cd $out

# Merge the sorted alignment files
echo "Merging sorted aligment file"
samtools merge -o merged.coordSorted.bam *.coordSorted.bam --threads $SLURM_CPUS_PER_TASK

# Index the coordinate-sorted bam file
echo "Indexing the merged bam"
samtools index merged.coordSorted.bam --threads $SLURM_CPUS_PER_TASK

echo "Done"
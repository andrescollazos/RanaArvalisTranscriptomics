#!/bin/bash
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Eval_Read_Rep_Transcriptome
#SBATCH -t 6:00:00
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
# bowtie2-build $transcriptome Trinity.fasta --threads $SLURM_CPUS_PER_TASK

# Perform the alignment to capture the read alignment statistics.
samples="$DIR/data/trinity/samples.txt"
while IFS=$'\t' read -r condition sample_id left right || [[ -n "$condition" ]]; do
    echo "------------------------------------------------------------------------"
    echo "Aligning ${sample_id}:"
    echo "------------------------------------------------------------------------"
    bowtie2 --threads $SLURM_CPUS_PER_TASK -q --no-unal -k 20 \
    -x Trinity.fasta \
    -1 "$DIR/data/trimmed_reads/$sample_id/$left" \
    -2 "$DIR/data/trimmed_reads/$sample_id/$right"  \
    2>"${sample_id}_align_stats.txt" | samtools view --threads $SLURM_CPUS_PER_TASK -Sb -o "${sample_id}.bam"

done < "$samples"

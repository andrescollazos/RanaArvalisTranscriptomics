#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Blast_Chunk
#SBATCH -t 4:00:00
#SBATCH -p node
#SBATCH --array=0-48%25
#SBATCH --cpus-per-task=20
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools blast/2.16.0+
source ../../.env

dir="$DIR/results/03_evaluation/02_count_full_length_transcripts"
cd $dir

sequences=(chunks/*.tmp.fasta)
query_file=${sequences[$SLURM_ARRAY_TASK_ID]}

echo "Performing $SLURM_ARRAY_TASK_ID BLASTx search"
blastx -query "$query_file" \
    -db swissprot \
    -out "blast.${SLURM_ARRAY_TASK_ID}.outfmt6" \
    -evalue 1e-20 \
    -outfmt 6 \
    -max_target_seqs 1 \
    -num_threads $SLURM_CPUS_PER_TASK
echo "Done"

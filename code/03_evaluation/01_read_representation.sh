source ../../.env
samples="$DIR/data/trinity/samples.txt"
out="$DIR/results/03_evaluation/00_read_representation"

# Dispatch an evaluation job for each
while IFS=$'\t' read -r condition sample_id left right || [[ -n "$condition" ]]; do
    echo "------------------------------------------------------------------------"
    echo "Dispatching alignment for ${sample_id}:"
    echo "------------------------------------------------------------------------"
    (echo   "#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J \"Align_$sample_id\"
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 2:30:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools bowtie2/2.5.2 samtools/1.20
source ../../.env

cd $out
#Perform the alignment to capture the read alignment statistics.
bowtie2 --threads \$SLURM_CPUS_PER_TASK -q --no-unal -k 20 \
    -x Trinity.fasta \
    -1 \"$DIR/data/trimmed_reads/$sample_id/$left\" \
    -2 \"$DIR/data/trimmed_reads/$sample_id/$right\"  \
    2>\"${sample_id}_align_stats.txt\" | samtools view --threads \$SLURM_CPUS_PER_TASK -Sb -o \"${sample_id}.bam\"

exit 0"
    ) | sbatch
    echo "${sample_id} Dispatched successfully"

done < "$samples"
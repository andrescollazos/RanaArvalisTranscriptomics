#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Prep_IGV_Launcher
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 2:30:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools samtools/1.20
source ../../.env

samples="$DIR/data/trinity/samples.txt"
out="$DIR/results/03_evaluation/00_read_representation"

cd $out
# Index the Trinity.fasta file
samtools faidx Trinity.fasta

cd -

# Dispatch an evaluation job for each sample
while IFS=$'\t' read -r condition sample_id left right || [[ -n "$condition" ]]; do
    echo "------------------------------------------------------------------------"
    echo "Dispatching IGV preparation for ${sample_id}:"
    echo "------------------------------------------------------------------------"
    (echo   "#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J \"Sort_${sample_id}_BAM\"
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

cd $out
# Sort the alignment by coordinate
samtools sort \"${sample_id}.bam\" -o \"${sample_id}.coordSorted.bam\" -@ \$SLURM_CPUS_PER_TASK

exit 0") | sbatch
    echo "${sample_id} Dispatched successfully"s

done < "$samples"
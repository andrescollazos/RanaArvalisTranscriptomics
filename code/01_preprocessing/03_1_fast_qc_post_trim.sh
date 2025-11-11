source ../../.env

out="$DIR/results/01_preprocessing/03_qc_trimmed_reads"
mkdir -p $out
cores=16

(echo "#!/bin/bash
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J FastQC_trimmed_reads
#SBATCH -e \"$out/%x.%j.er\"
#SBATCH -o \"$out/slurm-%x.%j.out\"
#SBATCH -t 3:00:00
#SBATCH -n $cores
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load FastQC/0.11.9

reads=\"$DIR/data/trimmed_reads/P32262_*/*.fastq.gz\"

fastqc \$reads -o $out -t $cores
exit 0") | sbatch
echo "QC analysis dispatched successfully"

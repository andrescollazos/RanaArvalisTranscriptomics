source ../../.env

fastqc_results="$DIR/results/01_preprocessing/03_qc_trimmed_reads"
out="$fastqc_results/multiqc"
mkdir -p $out
cores=2

(echo "#!/bin/bash
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J MultiQC
#SBATCH -e \"$out/%x.%j.er\"
#SBATCH -o \"$out/slurm-%x.%j.out\"
#SBATCH -t 30:00
#SBATCH -n $cores
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load MultiQC/1.22.2

multiqc $fastqc_results --outdir $out
exit 0") | sbatch
echo "MultiQC dispatched successfully"
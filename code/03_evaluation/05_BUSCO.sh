#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J BUSCO_Transcriptome
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 16:00:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools BUSCO/5.7.1
source ../../.env

out_dir="$DIR/results/03_evaluation/03_BUSCO"

cp "$DIR/data/transcriptome/Trinity.fasta" $SNIC_TMP
cd $SNIC_TMP

# Set up AUGUSTUS config locally
source $AUGUSTUS_CONFIG_COPY

# Run BUSCO
busco -i Trinity.fasta \
    -l $BUSCO_LINEAGE_SETS/tetrapoda_odb10/ \
    -o busco_out \
    -m transcriptome \
    -c $SLURM_CPUS_PER_TASK \
    -f

# Copy the output 
cp -r busco_out $out_dir
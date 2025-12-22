#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J TD2_LongOrfs
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 00:05:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load conda
export CONDA_ENVS_PATH=$PWD
conda activate TD2_env

source ../../.env
out="$DIR/results/07_annotation"
cd $out

TD2.LongOrfs \
  -t Trinity.fasta \
  --gene-trans-map Trinity.fasta.gene_trans_map \
  --threads ${SLURM_CPUS_PER_TASK}

conda deactivate

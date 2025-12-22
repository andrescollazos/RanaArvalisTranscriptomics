#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J TD2_Predict
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 16:00:00
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

TD2.Predict \
  -t Trinity.fasta \
  --retain-mmseqs-hits alnRes.m8

conda deactivate
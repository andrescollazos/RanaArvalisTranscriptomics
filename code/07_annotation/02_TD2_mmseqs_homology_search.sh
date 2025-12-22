#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J TD2_HomologySearch
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 00:20:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools MMseqs2/15-6f452
source ../../.env

dir="$DIR/results/07_annotation"
cd $dir

mmseqs easy-search \
  Trinity/longest_orfs.pep \
  $MMSEQS2_DATA/UniProtKB_Swiss-Prot \
  alnRes.m8 \
  tmp_mmseqs \
  -s 7.0 \
  --threads $SLURM_CPUS_PER_TASK
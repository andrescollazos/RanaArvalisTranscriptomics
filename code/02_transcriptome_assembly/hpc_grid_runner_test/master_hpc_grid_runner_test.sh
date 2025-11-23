#!/bin/bash
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J HPC_Grid_Runner_Master_Test
#SBATCH -t 3:00:00
#SBATCH -p core
#SBATCH --cpus-per-task=1
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

source ../../../.env

$HPC_GRID_RUNNER/hpc_cmds_GridRunner.pl -c ./test_commands.txt --grid_conf ./SLURM_test.conf


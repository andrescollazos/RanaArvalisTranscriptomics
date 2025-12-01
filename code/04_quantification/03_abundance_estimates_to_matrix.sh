#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J AlignEstAbundance
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 10:00:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

out="$DIR/results/04_quantification"
cd $out

# Get a list of the salmon quant.sf files
find . -maxdepth 2 -name "quant.sf" | tee salmon.quant_files.txt

singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/util/abundance_estimates_to_matrix.pl \
		--est_method salmon \
		--gene_trans_map Trinity.fasta.gene_trans_map \
        --quant_files salmon.quant_files.txt \
        --name_sample_by_basedir
#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J 05_Ex90N50
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 0:02:00
#SBATCH -p core
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

out="$DIR/results/04_quantification"
cd $out

echo "Compute Ex90N50"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/util/misc/contig_ExN50_statistic.pl \
        salmon.isoform.TMM.EXPR.matrix Trinity.fasta transcript | tee ExN50.transcript.stats

echo "Plot ExN50 statistic"
singularity exec --cleanenv \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg ls \
    /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript ExN50.transcript.stats
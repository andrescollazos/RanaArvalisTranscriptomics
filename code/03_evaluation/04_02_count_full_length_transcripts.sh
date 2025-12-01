#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Blast_Transcripts_Cov
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 4:00:00
#SBATCH -p node
#SBATCH --cpus-per-task=20
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools blast/2.16.0+
source ../../.env

dir="$DIR/results/03_evaluation/02_count_full_length_transcripts"
cd $dir

echo "Extracting data from the swissprot DB"
blastdbcmd \
  -db swissprot \
  -dbtype prot \
  -entry all \
  -out $dir/swissprot.fasta
echo "Done"

# Examine the percent of the target being aligned to by the best matching Trinity transcript
echo "Running analyze_blastPlus_topHit_coverage"
singularity exec --cleanenv \
    --env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	/usr/local/bin/util/analyze_blastPlus_topHit_coverage.pl blastx.first.outfmt6 Trinity.fasta swissprot.fasta
echo "Done"


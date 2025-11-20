#!/bin/bash
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Trinity_Normalization
#SBATCH -t 96:00:00
#SBATCH -p node
#SBATCH --mem=120G
#SBATCH --cpus-per-task=20
#SBATCH --error=%x.er
#SBATCH --output=slurm-%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
source ../../.env

# Copy the reads to the scratch
samples="$DIR/data/trinity/samples.txt"
while IFS=$'\t' read -r condition sample_id left right || [[ -n "$condition" ]]; do
    cp "$DIR/data/trimmed_reads/$sample_id/$left" "$SNIC_TMP/$left"
    cp "$DIR/data/trimmed_reads/$sample_id/$right" "$SNIC_TMP/$right"
done < "$samples"
cp $samples $SNIC_TMP/samples.txt

cd $SNIC_TMP
singularity exec --cleanenv	\
	--env LANG=C	\
	--env LC_ALL=C	\
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg	\
	Trinity --seqType fq	\
	--SS_lib_type FR \
	--max_memory 120G	\
	--samples_file samples.txt	\
	--output trinity_out	\
	--CPU 20	\
	--no_run_inchworm

# copy output in current directory
cp -Pr trinity_out/ $DIR/results/02_transcriptome_assembly

# TRINITY_SINGULARITY points to the path where Trinity v2.15.2 singularity image is located
# Trinity v2.15.2 singularity image was fetched from https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/
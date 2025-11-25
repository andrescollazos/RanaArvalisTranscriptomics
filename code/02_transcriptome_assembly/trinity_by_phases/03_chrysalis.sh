#!/bin/bash
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Trinity_Chrysalis
#SBATCH -t 5-00:00:00
#SBATCH -p node
#SBATCH --mem=120G
#SBATCH --cpus-per-task=20
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
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

# Copy the in silico normalization out
cp -rP "$DIR/results/02_transcriptome_assembly/trinity_out" $SNIC_TMP

cd $SNIC_TMP
tree -sh > "$DIR/results/02_transcriptome_assembly/pre_chrysalis_files.txt"
singularity exec --cleanenv \
	--env LANG=C \
	--env LC_ALL=C \
	$TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg \
	Trinity --seqType fq \
	--SS_lib_type FR \
	--max_memory 120G \
	--samples_file samples.txt \
	--output trinity_out \
	--CPU 20 \
	--no_distributed_trinity_exec

# Fix absolute paths in the commands for the Butterfly step
# Change internal path to Trinity to the actual path
sed -i 's|/usr/local/bin/util/support_scripts/../../Trinity|singularity exec --cleanenv --env LANG=C --env LC_ALL=C '"$TRINITY_SINGULARITY"'/trinityrnaseq.v2.15.2.simg Trinity|g' trinity_out/recursive_trinity.cmds
# Change path to the partition reads to the SCRATCH_PATH flag to be replaced in Butterfly
sed -i "s|$SNIC_TMP|SCRATCH_PATH|g" trinity_out/recursive_trinity.cmds

# Copy the output
cp -Pr trinity_out/ $DIR/results/02_transcriptome_assembly

# TRINITY_SINGULARITY points to the path where Trinity v2.15.2 singularity image is located
# Trinity v2.15.2 singularity image was fetched from https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/
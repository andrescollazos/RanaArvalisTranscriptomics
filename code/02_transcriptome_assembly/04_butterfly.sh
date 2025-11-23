#!/bin/bash
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Trinity_Butterfly
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
echo "Copying the samples..."
samples="$DIR/data/trinity/samples.txt"
while IFS=$'\t' read -r condition sample_id left right || [[ -n "$condition" ]]; do
	cp "$DIR/data/trimmed_reads/$sample_id/$left" "$SNIC_TMP/$left"
	cp "$DIR/data/trimmed_reads/$sample_id/$right" "$SNIC_TMP/$right"
done < "$samples"
cp $samples $SNIC_TMP/samples.txt

# Copy the in chrysalis out
echo "Copying the chrysalis out..."
cp -rP "$DIR/results/02_transcriptome_assembly/trinity_out" $SNIC_TMP

cd $SNIC_TMP
# Set the correct path to the partition reads in the scratch
sed -i "s|SCRATCH_PATH|$SNIC_TMP|g" trinity_out/recursive_trinity.cmds

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
	--grid_exec "$HPC_GRID_RUNNER/hpc_cmds_GridRunner.pl --grid_conf $SLURM_SUBMIT_DIR/SLURM.conf -c" \
	--full_cleanup \

tree -sh > "$DIR/results/02_transcriptome_assembly/final_files.txt"
# Copy the output
cp -Pr trinity_out/ $DIR/results/02_transcriptome_assembly

# TRINITY_SINGULARITY points to the path where Trinity v2.15.2 singularity image is located
# Trinity v2.15.2 singularity image was fetched from https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/
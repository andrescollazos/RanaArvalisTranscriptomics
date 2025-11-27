#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J Chunk_Transcripts
#SBATCH -e %x.%j.er
#SBATCH -o %x.%j.out
#SBATCH -t 0:20:00
#SBATCH -p core
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --error=%x.%j.er
#SBATCH --output=%x.%j.out
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

source ../../.env

dir="$DIR/results/03_evaluation/02_count_full_length_transcripts"
cd $dir

# If blast was previously executed
# cut -f1 blastx.outfmt6 | sort -u > done_sequences.txt


mkdir -p chunks
awk -v max=10000 '   
    BEGIN {
        # load all processed IDs into a hash
        while ((getline < "done_sequences.txt") > 0) done[$1] = 1
    }

    /^>/ {
        # extract Trinity ID: everything after ">" until first space
        split($1, a, " ")
        header = substr(a[1], 2)   # remove leading ">"

        # skip if this transcript already has BLAST results
        if (header in done) { skip = 1; next }
        skip = 0

        # create new chunk file every max transcripts
        if (count % max == 0) {
            if (out) close(out)
            out = sprintf("chunks/chunk_%03d.tmp.fasta", ++n)
        }
        count++
    }

    # print only if not skipped
    !skip { print > out }

' Trinity.fasta

echo "You must run 04_01_blast_transcripts_array.sh with --array=0-$(( $(ls chunks/*.tmp.fasta | wc -l) - 1 ))%20"
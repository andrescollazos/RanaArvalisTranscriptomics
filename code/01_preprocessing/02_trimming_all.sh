source ../../.env

base_in="${DIR}/data/raw_reads"
base_out="${DIR}/data/trimmed_reads"

 # First run: test script and measure time
# for d in "$base_in"/P32262_201; do
#    sample=$(basename "$d")

# All reads
# for d in "$base_in"/P32262_*; do
#    sample=$(basename "$d")

# Specific set of samples defined in samples.txt
for sample in $(cat samples.txt); do 
    d="${base_in}/${sample}"

    out_sample="$base_out/$sample"
    mkdir -p $out_sample
    cores=8
    input_R1=$(ls "$d"/*_R1_*.fastq.gz)
    input_R2=$(ls "$d"/*_R2_*.fastq.gz)
    
    # Skip if trimmed files already exist
    if ls "$out_sample"/*.fastq.gz >/dev/null 2>&1; then
        echo "Skipping ${sample} - trimmed files already exist."
        continue
    fi

    echo echo "Dispatching: ${sample}"
    (echo   "#!/bin/bash -l
#SBATCH -A uppmax2025-2-221
#SBATCH -M rackham
#SBATCH -J \"trim_$sample\"
#SBATCH -e \"$out_sample/%x.%j.er\"
#SBATCH -o \"$out_sample/slurm-%x.%j.out\"
#SBATCH -t 45:00
#SBATCH -n $cores
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load trimmomatic/0.39

input_R1=\"$input_R1\"
input_R2=\"$input_R2\"
out_R1_paired=\"$out_sample/${sample}_R1_paired.fastq.gz\"
out_R1_upaired=\"$out_sample/${sample}_R1_upaired.fastq.gz\"
out_R2_paired=\"$out_sample/${sample}_R2_paired.fastq.gz\"
out_R2_upaired=\"$out_sample/${sample}_R2_upaired.fastq.gz\"

trimmomatic PE -threads $cores -phred33 \\
    \$input_R1 \$input_R2 \\
    \$out_R1_paired \$out_R1_upaired \$out_R2_paired \$out_R2_upaired \\
    ILLUMINACLIP:\$TRIMMOMATIC_ROOT/adapters/TruSeq3-PE.fa:2:30:10 \\
    SLIDINGWINDOW:4:20 \\
    TRAILING:3 \\
    MINLEN:50
exit 0"
    ) | sbatch
    echo "${sample} Dispatched successfully"

done;

# $TRIMMOMATIC_ROOT is automatically loaded by "module load trimmomatic/0.39"
source ../../.env
dir="$DIR/results/03_evaluation/00_read_representation"
cd $dir

echo -e "Sample\tTotal reads\tOverall alignment (%)\tConcordant = 1 (%)\tConcordant > 1 (%)\tOther Categories (%)\tUnmapped (%)" > align_stats_summary.txt

for f in *_align_stats.txt; do
    sample=$(basename "$f" _align_stats.txt)

    total_reads=$(grep "reads; of these" "$f" | awk '{print $1}')
    overall_align=$(grep -Eo "[0-9]+\.[0-9]+% overall alignment rate" "$f" | cut -d% -f1)

    concord1=$(grep "aligned concordantly exactly 1 time" "$f" | awk -F'[()%]' '{print $2}')
    concordGT1=$(grep "aligned concordantly >1 times" "$f" | awk -F'[()%]' '{print $2}')

    unmapped=$(echo "100 - $overall_align" | bc)
    others=$(echo "$overall_align - $concord1 - $concordGT1" | bc)

    echo -e "${sample}\t${total_reads}\t${overall_align}\t${concord1}\t${concordGT1}\t${others}\t${unmapped}" >> align_stats_summary.txt
done
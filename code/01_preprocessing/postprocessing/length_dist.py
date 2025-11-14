import sys
import csv

# $RESULTS/01_preprocessing/03_qc_trimmed_reads/comparison/length_distribution.txt
path = sys.argv[1]

values = []

with open(path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        ratio = int(row["reads_longer_51"]) / int(row["total_reads"])
        values.append(ratio)

avg = sum(values) / len(values)
print("Avg:", avg * 100)
import sys
import json

# $RESULTS/01_preprocessing/03_qc_trimmed_reads/comparison/counts.json
path = sys.argv[1]

# load JSON
with open(path) as f:
    data = json.load(f)

raw = data["raw"]
trimmed = data["trimmed"]

diff = {}

for sample in raw:
    unique_var = 1 - (trimmed[sample][0] / raw[sample][0])
    raw_duplicate_level = raw[sample][1] / (raw[sample][1] + raw[sample][0])
    trimmed_duplicate_level = trimmed[sample][1] / (trimmed[sample][1] + trimmed[sample][0])
    duplicate_var = trimmed_duplicate_level - raw_duplicate_level

    diff[sample] = [unique_var, duplicate_var]

# print(diff)

# average unique reduction
total_unique_reduction = [v[0] for v in diff.values()]
unique_avg_reduction = sum(total_unique_reduction) / len(total_unique_reduction)
print("Average Unique Reduction:", f"{unique_avg_reduction * 100:.4f} %")

# average duplicate diff
total_duplicate_reduction = [v[1] for v in diff.values()]
duplicate_avg_diff = sum(total_duplicate_reduction) / len(total_duplicate_reduction)
print("Average Duplicate Diff:", f"{duplicate_avg_diff * 100:.4f} %")
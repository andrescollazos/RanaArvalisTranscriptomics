import sys
import json

# $RESULTS/01_preprocessing/03_qc_trimmed_reads/comparison/overrepresentation.json
path = sys.argv[1]

# load JSON
with open(path) as f:
    data = json.load(f)

raw = data["raw"]
trimmed = data["trimmed"]

# compute diff
diff = {
    sample: [
        trimmed[sample][0] - raw[sample][0],
        trimmed[sample][1] - raw[sample][1]
    ]
    for sample in raw
}

# print(diff)

# average top
top_seq = [v[0] for v in diff.values()]
top_avg = sum(top_seq) / len(top_seq)
print("Average top:", f"{top_avg * 100} %")

# average remaining
remain_seq = [v[1] for v in diff.values()]
remain_avg = sum(remain_seq) / len(remain_seq)
print("Average remaining:", f"{remain_avg * 100} %")
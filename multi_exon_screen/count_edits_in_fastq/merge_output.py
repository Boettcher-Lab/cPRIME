import pandas as pd
import sys

# function definitions
def read_file(path, sep = "\t"):
    result = {}
    with open(path, "r") as f:
        for line in f:
            if line.strip() == "" or line.startswith('#'):
                continue
            line = line.strip()
            parts = line.split(sep)

            count = int(parts[0])
            search_sequence = parts[1]

            result[search_sequence] = count

    return result

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


# argument parser
forward_count_path = sys.argv[1]
reverse_count_path = sys.argv[2]
ss_path = sys.argv[3]
out_path = sys.argv[4]

search_seqs = pd.read_csv(ss_path, sep=r'\s+', comment = '#', header = None)
search_seqs.columns = ["ss_id", "seq"]
search_seqs = search_seqs["seq"]


forward = pd.read_csv(forward_count_path, sep = "\t")
assert len(forward["total_reads"].unique()) == 1
total_reads = forward["total_reads"][0]

reverse = None
if reverse_count_path != "NONE":
    reverse = pd.read_csv(reverse_count_path, sep = "\t")
    reverse["seq"] = [reverse_complement(seq) for seq in reverse["seq"]]
    assert len(reverse["total_reads"].unique()) == 1

if reverse is not None:
    assert reverse["total_reads"].unique() == forward["total_reads"].unique()
    merged = pd.merge(forward[["seq", "counts"]], reverse[["seq", "counts"]], on='seq', how='outer', suffixes=["_fwd", "_rev"])
    merged["counts_complete"] = merged["counts_fwd"] + merged["counts_rev"]
    merged["total_reads"] = total_reads
else:
    merged = forward.copy()
    merged.columns = ["seq", "counts_complete", "total_reads"]

merged = merged.merge(search_seqs, on = "seq", how = "outer")
merged = merged.fillna(0)
merged = merged.sort_values("counts_complete", ascending=False)
merged["fraction"] = (merged["counts_complete"] + 1) / (merged["total_reads"] + 1)
merged["RPM"] = merged["fraction"] * 1000000

merged.to_csv(out_path, sep = '\t', index = False)
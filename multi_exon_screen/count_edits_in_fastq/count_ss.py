import sys
import gzip
from collections import OrderedDict
import ahocorasick



def load_motifs(motif_path, motif_col=2, sep=None):
    motifs = []

    with open(motif_path, "rt") as motif_file:
        for line in motif_file:
            line = line.strip("\n")
            if not line.strip():
                continue

            fields = line.split(sep)
            if len(fields) < motif_col:
                raise ValueError(f"not enough columns in line: {line}")

            motif = fields[motif_col - 1].strip().upper()
            if not motif:
                sys.stderr.write(
                    f"WARNING: skipping empty motif in line: {line}\n"
                )
                continue

            motifs.append(motif)

    return list(set(motifs)) # make unique


def build_automaton(motifs):
    A = ahocorasick.Automaton()
    for i, motif in enumerate(motifs):
        A.add_word(motif, (i, motif))
    A.make_automaton()
    return A


def count_motifs_in_fastq(fastq_file, motifs, automaton):
    counts = {m: 0 for m in motifs}
    total_reads = 0

    with gzip.open(fastq_file, "rt") as infile:
        line_num = 0
        for line in infile:
            line_num += 1
            if line_num % 4 != 2: # skip non sequence lines
                continue

            seq = line.strip().upper()
            total_reads += 1

            motifs_in_read = set()
            for _, (idx, motif) in automaton.iter(seq):
                motifs_in_read.add(motif)
            for motif in motifs_in_read:
                counts[motif] += 1

    return total_reads, counts



if len(sys.argv) != 4:
    sys.stderr.write(
        "Need three arguments: motif_file_path, fastq_input_path and output_path. The motif_file should be in the format: motif_id <space> motif"
    )
    sys.exit(1)

motif_file_path = sys.argv[1]
fastq_input_path = sys.argv[2]
output_path = sys.argv[3]

# load the motifs
motifs = load_motifs(motif_file_path, motif_col=2)
if not motifs:
    raise ValueError("no motifs found. please check input motif file.")

# count motifs in file using the ahocorasick algorithm
automaton = build_automaton(motifs)
total_reads, counts = count_motifs_in_fastq(fastq_input_path, motifs, automaton)

# write output
with open(output_path, "wt") as out:
    out.write("seq\tcounts\ttotal_reads\n")
    for motif in motifs:
        c = counts[motif]
        out.write('\t'.join([str(motif), str(c), str(total_reads)]) + "\n")


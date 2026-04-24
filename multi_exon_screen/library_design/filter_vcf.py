import sys
import argparse
import gzip


parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", "--input",  default="", help="path to input.vcf file")
parser.add_argument("-o", "--output", default="", help="output file path. If not given will default to stdout")
parser.add_argument("-c", "--chromosome", help="The chromosome of interest. Must be an exact match (provide chr prefix if needed)")
parser.add_argument("-s", "--start", help="The start position of the region of interest")
parser.add_argument("-e", "--end", help="The end position of the region of interest")

args = parser.parse_args()

if args.output != "":
    sys.stdout = open(args.output, 'w')

input_path = args.input
if input_path != "":
    if input_path.endswith(".gz"):
        input_file = gzip.open(input_path, "rt")
    else:
        input_file = open(input_path, 'r')
else:
    input_file = sys.stdin

chrom_oi = args.chromosome
start_oi = int(args.start)
end_oi = int(args.end)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

for line in input_file:

    if line.strip() == '': # skip empty lines
        continue

    if line.startswith('#'): # preserve comments and header lines
        print(line.strip('\n'))
        continue

    parts = line.split('\t') #CHROM POS      ID         REF   ALT    QUAL  FILTER   INFO

    chrom = parts[0]
    pos = int(parts[1])

    if chrom == chrom_oi and start_oi <= pos <= end_oi:
        print(line.strip('\n'))

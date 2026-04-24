# python3 script_name.py input_fastq.gz output_masked_fastq.gz
import sys
import gzip

# Define the function to replace bases
def replace_bases(seq, qual):
    new_seq = ''
    for i in range(len(seq)):
        if qual[i] in ['F', 'G', 'H', 'I', 'J', 'K']:
            new_seq += seq[i]
        else:
            new_seq += 'X'
    return new_seq

# Read in the input fastq file and replace bases
with gzip.open(sys.argv[1], 'rt') as f_in, gzip.open(sys.argv[2], 'wt') as f_out:
    while True:
        # Read in the four lines of each record
        header = f_in.readline().strip()
        if not header:
            break  # end of file
        seq = f_in.readline().strip()
        plus = f_in.readline().strip()
        qual = f_in.readline().strip()

        # Replace the bases in the sequence
        masked_seq = replace_bases(seq, qual)

        # Write the masked record to the output file
        f_out.write(header + '\n')
        f_out.write(masked_seq + '\n')
        f_out.write(plus + '\n')
        f_out.write(qual + '\n')

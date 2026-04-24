import sys
import gzip
import numpy as np

# constants
QMIN = ord('F')
MASK = ord('X')
NEWLINE = ord('\n')
CHUNK_SIZE = 64 * 1024 * 1024


# functions
def process_fastq_chunk(buf: bytes) -> bytes:
    arr = np.frombuffer(buf, dtype=np.uint8)

    nl = np.flatnonzero(arr == NEWLINE)
    n_lines = len(nl)

    if n_lines % 4 != 0:
        raise ValueError("number of lines in chunk is not dividable by 4. Please check input FASTQ.")

    starts = np.empty(n_lines, dtype=np.int64)
    starts[0] = 0
    starts[1:] = nl[:-1] + 1
    ends = nl

    seq_starts = starts[1::4]
    seq_ends = ends[1::4]
    qual_starts = starts[3::4]
    qual_ends = ends[3::4]

    out = np.frombuffer(bytearray(buf), dtype=np.uint8)

    for seq_start, seq_end, qual_start, qual_end in zip(seq_starts, seq_ends, qual_starts, qual_ends):
        if (seq_end - seq_start) != (qual_end - qual_start):
            raise ValueError(
                f"length of quality and sequence strings are unequal: seq={seq_end-seq_start} while qual={qual_end-qual_start}"
            )

        q = arr[qual_start:qual_end]
        bad = q < QMIN
        out[seq_start:seq_end][bad] = MASK

    return out.tobytes()

# start script

if len(sys.argv) != 3:
    sys.stderr.write(
        "Need two arguments: in_path and out_path"
    )
    sys.exit(1)

in_path = sys.argv[1]
out_path = sys.argv[2]
outfile = gzip.open(out_path, "wb")

carry = b""
with gzip.open(in_path, "rb") as fin:
    while True:
        block = fin.read(CHUNK_SIZE)
        if not block:
            break

        data = carry + block
        arr = np.frombuffer(data, dtype=np.uint8)
        nl = np.flatnonzero(arr == NEWLINE)
            
        complete_lines = (len(nl) // 4) * 4
        if complete_lines == 0:
            carry = data
            continue

        cut = nl[complete_lines - 1] + 1
        chunk = data[:cut]
        carry = data[cut:]

        outfile.write(process_fastq_chunk(chunk))

    # handle remaining bytes
    if carry:
        arr = np.frombuffer(carry, dtype=np.uint8)
        nl = np.flatnonzero(arr == NEWLINE)

        if len(nl) == 0:
            raise ValueError("input file ended with incomplete fastq record")

        if len(nl) % 4 != 0 or nl[-1] != len(carry) - 1:
            raise ValueError("input file ended with incomplete fastq record")

        outfile.write(process_fastq_chunk(carry))

outfile.close()

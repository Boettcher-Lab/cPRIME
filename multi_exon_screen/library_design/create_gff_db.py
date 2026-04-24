import gffutils
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('input', help = 'Ensembl gff file')
ap.add_argument('database', help = 'Database to be created')
ap.add_argument('--testN', type = int, help = 'run a test using only the first N features, and then print out some example feature IDs and their attributes')
ap.add_argument('--force', action = 'store_true', help = 'Overwrite existing database')
args = ap.parse_args()


def first_n_features(data, n=5000):
    for i, feature in enumerate(gffutils.iterators.DataIterator(data)):
        if i > n:
            break
        yield feature


id_spec = {
    'exon': 'exon_id',
    'gene': 'gene_id',
    'transcript': 'transcript_id',
    'mRNA': 'transcript_id',
}

if args.testN is None:
    data = args.input
else:
    data = first_n_features(args.input, args.testN)

db = gffutils.create_db(
    data,
    args.database,
    id_spec=id_spec,
    disable_infer_genes=True,
    disable_infer_transcripts=True,
    merge_strategy='create_unique',
    verbose=True,
    force=args.force,
)

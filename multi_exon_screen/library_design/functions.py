import re

def translate_codon(codon):
    codon = codon.upper().replace("T", "U")
    codon2aa={'UUU':'Phe',
			'UUC':'Phe',
			'UUA':'Leu',
			'UUG':'Leu',
			'CUU':'Leu',
			'CUC':'Leu',
			'CUA':'Leu',
			'CUG':'Leu',
			'AUU':'Ile',
			'AUC':'Ile',
			'AUA':'Ile',
			'AUG':'Met',
			'GUU':'Val',
			'GUC':'Val',
			'GUA':'Val',
			'GUG':'Val',
			'UCU':'Ser',
			'UCC':'Ser',
			'UCA':'Ser',
			'UCG':'Ser',
			'CCU':'Pro',
			'CCC':'Pro',
			'CCA':'Pro',
			'CCG':'Pro',
			'ACU':'Thr',
			'ACC':'Thr',
			'ACA':'Thr',
			'ACG':'Thr',
			'GCU':'Ala',
			'GCC':'Ala',
			'GCA':'Ala',
			'GCG':'Ala',
			'UAU':'Tyr',
			'UAC':'Tyr',
			'UAA':'ter',
			'UAG':'ter',
			'CAU':'His',
			'CAC':'His',
			'CAA':'Gln',
			'CAG':'Gln',
			'AAU':'Asn',
			'AAC':'Asn',
			'AAA':'Lys',
			'AAG':'Lys',
			'GAU':'Asp',
			'GAC':'Asp',
			'GAA':'Glu',
			'GAG':'Glu',
			'UGU':'Cys',
			'UGC':'Cys',
			'UGA':'ter',
			'UGG':'Trp',
			'CGU':'Arg',
			'CGC':'Arg',
			'CGA':'Arg',
			'CGG':'Arg',
			'AGU':'Ser',
			'AGC':'Ser',
			'AGA':'Arg',
			'AGG':'Arg',
			'GGU':'Gly',
			'GGC':'Gly',
			'GGA':'Gly',
			'GGG':'Gly'}
    return codon2aa[codon].lower()


def get_codons(amino_acid):
    aa2codon = {
        "ala":["GCU", "GCC", "GCA", "GCG"],
        "arg":["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "asn":["AAU", "AAC"],
        "asp":["GAU", "GAC"],
        "asx":[],
        "cys":["UGU", "UGC"],
        "glu":["GAA", "GAG"],
        "gln":["CAA", "CAG"],
        "glx":[],
        "gly":["GGU", "GGC", "GGA", "GGG"],
        "his":["CAU", "CAC"],
        "ile":["AUU", "AUC", "AUA"],
        "leu":["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
        "lys":["AAA", "AAG"],
        "met":["AUG"],
        "phe":["UUU", "UUC"],
        "pro":["CCU", "CCC", "CCA", "CCG"],
        "ser":["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
        "thr":["ACU", "ACC", "ACA", "ACG"],
        "trp":["UGG"],
        "tyr":["UAU", "UAC"],
        "val":["GUU", "GUC", "GUA", "GUG"],
        "ter":["UAA", "UAG", "UGA"]
    }
    return [x.replace('U', "T") for x in aa2codon[amino_acid]]


# functions to calculate synonymous mutations
def genomic2cdna(exons, genomic_position):
    var_exon_info = exons[(exons["start"] <= genomic_position) & (genomic_position <= exons["end"])] # extract the exon where the variant is located
    if len(var_exon_info) == 0:
        return None
    var_exon_info = var_exon_info.reset_index()
    var_dist_to_exon_start = genomic_position - var_exon_info.loc[0,"start"]
    var_cdna_pos = var_dist_to_exon_start + var_exon_info.loc[0,"cdna_start"]
    return var_cdna_pos


def cdna2codon(cdna_pos):
    # compute codon position
    # consider a codon:
    # A T G
    # ^ ^ ^
    # 0 1 2
    return (cdna_pos - 1) % 3


def compute_adjacent_codons(genom_pos, codon_pos, n_codons=1, verbose = False):
    dist_to_next_codon_start = 3 - codon_pos
    dist_to_prev_codon_start = codon_pos + 3

    if verbose:
        print(dist_to_next_codon_start)
        print(dist_to_prev_codon_start)

    next_codon_start = genom_pos + dist_to_next_codon_start
    prev_codon_start = genom_pos - dist_to_prev_codon_start

    result = [next_codon_start, prev_codon_start]

    old_next_codon_start = next_codon_start
    old_prev_codon_start = prev_codon_start
    for n in range(n_codons-1):
        old_next_codon_start += 3
        old_prev_codon_start -= 3
        result.extend([
            old_next_codon_start,
            old_prev_codon_start
        ])
    
    return result


def compute_codon_start(genom_pos, codon_pos):
    return genom_pos - codon_pos


def get_required_mutations(chrom, ref_seq, alt_seq, pos):
    mutations_req = []
    for i, alt_base in enumerate(alt_seq):
        ref_base = ref_seq[i]
        if ref_base != alt_base:
            mutations_req.append((chrom, pos + i, ref_base, alt_base))
    return mutations_req


def get_synonymous_variants(codon, chrom, genomic_codon_start):
    aa = translate_codon(codon)
    syn_codons = get_codons(aa)
    variants = [get_required_mutations(chrom, codon, syn_codon, genomic_codon_start) for syn_codon in syn_codons]
    variants = [x[0] for x in variants if len(x) == 1]
    return variants


# functions to compute the mutated sequence
def introduce_mutation(seq, pos, ref, alt, format = "primedesign"):
    if format == "primedesign":
        mutation = "(" + ref + "/" + alt + ")"
    elif format == "alt":
        mutation = alt
    offset = len(mutation) - len(ref)
    return seq[:pos -1] + mutation + seq[pos:], offset


def get_mutated_sequence(mutations, genome, format = "primedesign", verbose = False, flanking_dist = 300): # input: a list of lists with each element containing information about a variant in format: (chrom, pos, ref, alt)
    all_positions = [int(x[1]) for x in mutations]
    chrom = mutations[0][0] # assume they are all equal!
    
    wt_seq_start =  min(all_positions)
    wt_seq_end = max(all_positions)
    if format == "primedesign":
        wt_seq_start -= flanking_dist
        wt_seq_end += flanking_dist
    wt_seq = genome.fetch(chrom, wt_seq_start - 1, wt_seq_end)
    mut_seq = wt_seq

    offset = 0
    last_pos = None
    for mutation in mutations:
        pos = int(mutation[1])
        ref = mutation[2]
        alt = mutation[3]

        if last_pos is None or last_pos < pos:
            mut_seq, offset = introduce_mutation(mut_seq, (pos + offset) - wt_seq_start + 1, ref, alt, format = format)
        else:
            mut_seq, offset = introduce_mutation(mut_seq, pos - wt_seq_start + 1, ref, alt, format = format)
        last_pos = pos

        if verbose:
            print(mut_seq)
    return chrom, wt_seq_start, wt_seq, mut_seq


def find_between(s, prefix, postfix):
    res = re.search(prefix+r'(.*?)'+postfix, s)
    if res is not None:
        res = res.group(1)
    return res
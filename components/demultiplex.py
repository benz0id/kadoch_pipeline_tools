from typing import List

get_comp = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


def convert_to_rev_comp(seqs: List[str]) -> List[str]:
    """
    Converts each of seqs to its reverse complement.
    :param seqs: A list of DNA sequences.
    :return:
    """

    rev_comps = []
    for seq in seqs:
        comp = ''
        for nt in seq:
            comp += get_comp[nt]
        rev_comp = comp[::-1]
        rev_comps.append(rev_comp)
    return rev_comps


seqs = [
    "AAGAGGCA",
    "GTAGAGGA",
    "TGGATCTG",
    "CCGTTTGT",
    "TGCTGGGT",
    "AGGTTGGG",
    "GTGTGGTG",
    "TGGGTTTC"
]

print('\n'.join(convert_to_rev_comp(seqs)))
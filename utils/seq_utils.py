from typing import List

get_comp = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


def convert_to_rev_comps(seqs: List[str]) -> List[str]:
    """
    Converts each of seqs to its reverse complement.
    :param seqs: A list of DNA sequences.
    :return:
    """

    rev_comps = []
    for seq in seqs:
        rev_comps.append(get_rev_comp(seq))
    return rev_comps


def get_rev_comp(seq: str) -> str:
    """
    Convert the given seq to it's reverse compliment.
    :param seq: A series of non-degenerate nucleotides
    :return: The reverse compliment of <seq>.
    """
    comp = ''
    for nt in seq:
        comp += get_comp[nt]
    return comp[::-1]


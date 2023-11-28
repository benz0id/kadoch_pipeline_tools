import os
from typing import Union, List
from pathlib import Path

from utils.path_manager import cmdify


def fetch_ref_seqs(ids: Union[Path, List[str]], outfile: Path,
                   id_type: str = 'hgnc', db: str = 'protein',
                   organism: str = 'homo sapien') -> str:
    """
    Fetch reference sequences from the NCBI using the entrez system.
    :param file: Path to file containing some number of identifiers
        (one per line).
    :param id_type: The type of the given IDs.
    :param db: The type of sequences to fetch.
    :param outfile: File to write fastqs to.
    :param organism: The organism to select.
    """

    cmd = cmdify(f"esearch -db {db} -query",
                 f"\"{organism}[Organism] AND refseq[Filter] AND")

    if isinstance(ids, Path):
        with open(ids, 'r') as infile:
            for line in infile.readlines():
                if not line.strip():
                    continue
                cmd += ' ' + line.strip()

    elif isinstance(ids, list):
        for line in ids:
            if not line.strip():
                continue
            cmd += ' ' + line.strip()

    else:
        raise ValueError("Unrecognised input type")

    if id_type == 'hgnc':
        cmd += "[Gene Name]\""

    cmd += ' ' + cmdify('| efetch -format fasta',
                        '>', outfile)

    return cmd





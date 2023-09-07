import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from utils.fetch_files import get_matching_files
from utils.job_formatter import ExecParams
from utils.path_manager import cmdify

logger = logging.getLogger(__name__)


def merge_fastqs(fastq1_dir: Path, fastq2_dir: Path, mapping: Dict[str, str],
                 out_dir: Path, paired_end: bool = False,
                 verbose: bool = False) \
        -> Tuple[List[str], List[Path]]:
    """
        Merges the fastqs in <fastq1_dir> with the fastqs in <fastq2_dir>. Uses the
    <mapping> to search for files to merge.

        Combined fastqs are named using the name of the fastq1_dir fastq, and
    the information contained in the mapping.

        Does *not* create out directory.

    :param fastq1_dir: The directory containing the first set of fastq files.
    :param fastq2_dir: The directory contianing hte second set of fastq files.
    :param mapping: A dictionary that maps some identifiable substring of each
        fastq1 name to an identifiable substring of a fastq2 name.
    :param out_dir: The directory into which merged fastqs will be placed.
    :param paired_end: Whether the provided fastqs are paired end - will treat
        files with R1 and R2 in their names seperately.
    :param exec_params: Parameters with which to construct the merge job.
    :return: A list of commands that will create the merged fastqs once
        executed, and the paths to those resultant fastqs.


    Example:

        fastqs1
            1234_CG123_22_R1_001.fastq.gz
            1234_CG123_22_R2_001.fastq.gz
        fastqs2
            1234_CG500_22_R1_001.fastq.gz
            1234_CG500_22_R2_001.fastq.gz

    $ mkdir out_dir
    >>> cmds = merge_fastqs(Path("fastqs1"), Path("fastqs2"),
    >>>         {"CG123": "CG500"}, Path("out_dir"), paired_end=True)[0]
    >>> os.system(cmds[0])
    >>> os.system(cmds[1])

        out_dir
            1234_CG123&CG500_22_R1_001.fastq.gz
            1234_CG123&CG500_22_R2_001.fastq.gz

    """
    to_combine = []

    def get_fastq_pair(f1_patt: str, f2_patt: str,
                       f1_avoid: List[str] = None,
                       f2_avoid: List[str] = None) -> None:
        """
        Find a pair of fastq files and add them to <to_combine>.
        """

        # Get matching fastqs
        matching1 = get_matching_files(fastq1_dir, [f1_patt], f1_avoid,
                                       containing=True)
        matching2 = get_matching_files(fastq2_dir, [f2_patt], f2_avoid,
                                       containing=True)

        # Check that each identifier successfully found a single fastq.
        s = ''
        a1 = ''
        a2 = ''
        if f1_avoid:
            a1 = ' and avoiding ' + str(f1_avoid)
        if f2_avoid:
            a2 = ' and avoiding ' + str(f2_avoid)

        if len(matching1) == 0:
            s += f"Failed to find pattern matching \"{f1_patt}\"{a1} in {fastq1_dir}"
        if len(matching2) == 0:
            s += f"Failed to find pattern matching \"{f2_patt}\"{a2} in {fastq2_dir}"
        if len(matching1) > 1:
            mult = '\n\t'.join(matching1)
            s += f"Found multiple matches for \"{f1_patt}\"{a1} in {fastq1_dir}" \
                 f"\n{mult}"
        if len(matching2) > 1:
            mult = '\n\t'.join(matching2)
            s += f"Found multiple matches for \"{f2_patt}\"{a2} in {fastq2_dir}" \
                 f"\n{mult}"
        if s:
            print(s, sys.stderr)
            logger.fatal(s)
            raise ValueError(s)


        # Create the name for the new fastq file.
        fastq1_name = matching1[0]
        fastq2_name = matching2[0]
        name_comps = fastq1_name.split(f1_patt)

        new_name = name_comps[0] + f1_patt + '&' + f2_patt + name_comps[1]
        to_combine.append((out_dir / new_name,
                           fastq1_dir / fastq1_name,
                           fastq2_dir / fastq2_name))

    # Get all fastq pairs, with slightly different PE behaviour.
    if not paired_end:
        for id1 in mapping:
            id2 = mapping[id1]
            get_fastq_pair(id1, id2)
    else:
        for id1 in mapping:
            id2 = mapping[id1]
            get_fastq_pair(id1, id2, ['R2'], ['R2'])
            get_fastq_pair(id1, id2, ['R1'], ['R1'])

    # Convert all fastq pairs into cat commands.
    cmds = []
    new_fastqs = []
    cmd_str_rep = []

    for new_fastq, fastq1, fastq2 in to_combine:
        cmds.append(cmdify('cat', fastq1, fastq2, '>', new_fastq))
        cmd_str_rep.append(cmdify('cat', fastq1.name, fastq2.name,
                                  '>', new_fastq.name))
        new_fastqs.append(new_fastq)

    s = 'Creating the following merged fastqs.\n\t' + '\n\t'.join(cmd_str_rep)
    logger.info(s)
    if verbose:
        print(s)

    return cmds, new_fastqs



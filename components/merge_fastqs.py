import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from utils.fetch_files import get_matching_files, get_matching_strs
from utils.job_formatter import ExecParams
from utils.path_manager import cmdify
from utils.utils import ExperimentalDesign

logger = logging.getLogger(__name__)


def merge_fastqs_grouped(grouped_fastqs: List[List[Path]],
                         group_conditions: List[str], group_reps: List[str],
                         out_dir: Path) -> Tuple[List[str], List[Path]]:
    """
    For each group of fastq files in <grouped_fastqs>, combine all fastqs in
    that group into a single fastq file.

    :param grouped_fastqs: A list of grouped paths to fastq files. The fastqs
        in these groups will be combined into a single fastq.
    :param group_conditions: For each of <grouped_fastqs>, the condition under
        which all the given fastqs were conducted.
    :param group_reps: For each of the given <grouped_fastqs>, the replicate of
        that group.
    :param out_dir: The directory into which the combined fastqs will be
        placed.
    :return: A list of commands that would result in the creation of fastqs in
        the <out_dir> and a list of the resultant filepaths.
    """

    cmds = []
    paths = []

    if not len(grouped_fastqs) == len(group_conditions) == len(group_reps):
        raise ValueError(
            "grouped_fastqs, group_conditions, and group_reps must be of the same length. "
            "Currently:\n"
            f"grouped_fastqs: {len(grouped_fastqs)}\n"
            f"group_conditions: {len(group_conditions)}\n"
            f"group_reps: {len(group_reps)}")

    group_info = zip(grouped_fastqs, group_conditions, group_reps)

    for group, condition, rep in group_info:

        samples = []
        for fastq in group:
            samples.append(fastq.name.split('_')[1])

        tail = group[0].name.split('.')[0].split('_')[-2]
        if 'R1' in tail:
            read_str = '_R1'
        elif 'R2' in tail:
            read_str = '_R2'
        else:
            read_str = ''

        out_filename = '0MERGED0_' + '-'.join(samples) + '_' + \
                       condition + '_' + 'Rep' + str(rep) + read_str + \
                       '.fastq.gz'
        out_fastq = out_dir / out_filename

        cmds.append(cmdify('cat', *group, '>', out_dir / out_filename))
        paths.append(out_fastq)

    return cmds, paths


def merge_fastqs(fastqs1: List[Path], fastqs2: List[Path],
                 mapping: Dict[str, str], out_dir: Path,
                 paired_end: bool = False, verbose: bool = False) \
        -> Tuple[List[str], List[Path]]:
    """
        Merges the fastqs in <fastq1_dir> with the fastqs in <fastq2_dir>. Uses the
    <mapping> to search for files to merge.

        Combined fastqs are named using the name of the fastq1_dir fastq, and
    the information contained in the mapping.

        Does *not* create out directory.

    :param fastqs1: A list of fastq files.
    :param fastqs2: A seperate list of fastq files.
    :param mapping: A dictionary that maps some identifiable substring of each
        fastq1 name to an identifiable substring of a fastq2 name.
    :param out_dir: The directory into which merged fastqs will be placed.
    :param paired_end: Whether the provided fastqs are paired end - will treat
        files with R1 and R2 in their names seperately.
    :param exec_params: Parameters with which to construct the merge job.
    :return: A list of commands that will create the merged fastqs once
        executed, and the paths to those resultant fastqs.

    """
    to_combine = []

    def get_fastq_pair(f1_patt: str, f2_patt: str,
                       f1_avoid: List[str] = None,
                       f2_avoid: List[str] = None) -> None:
        """
        Find a pair of fastq files and add them to <to_combine>.
        """

        # Get matching fastqs
        matching1 = get_matching_strs([f.name for f in fastqs1], [f1_patt],
                                      f1_avoid, containing=True, inds=True)
        matching1 = [fastqs1[ind] for ind in matching1]
        matching2 = get_matching_strs([f.name for f in fastqs2], [f2_patt],
                                      f2_avoid, containing=True, inds=True)
        matching2 = [fastqs1[ind] for ind in matching2]

        # Check that each identifier successfully found a single fastq.
        s = ''
        a1 = ''
        a2 = ''
        if f1_avoid:
            a1 = ' and avoiding ' + str(f1_avoid)
        if f2_avoid:
            a2 = ' and avoiding ' + str(f2_avoid)

        if len(matching1) == 0:
            s += f"Failed to find pattern matching \"{f1_patt}\"{a1}"
        if len(matching2) == 0:
            s += f"Failed to find pattern matching \"{f2_patt}\"{a2}"
        if len(matching1) > 1:
            mult = '\n\t'.join([str(p) for p in matching1])
            s += f"Found multiple matches for \"{f1_patt}\"{a1}" \
                 f"\n{mult}"
        if len(matching2) > 1:
            mult = '\n\t'.join([str(p) for p in matching2])
            s += f"Found multiple matches for \"{f2_patt}\"{a2}" \
                 f"\n{mult}"
        if s:
            print(s, sys.stderr)
            logger.fatal(s)
            raise ValueError(s)

        # Create the name for the new fastq file.
        f1 = matching1[0]
        f2 = matching2[0]

        fastq1_name = f1.name

        name_comps = fastq1_name.split(f1_patt)

        new_name = '0MERGED0' + f1_patt + '-' + f2_patt + name_comps[1]
        to_combine.append((out_dir / new_name, f1, f2))

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

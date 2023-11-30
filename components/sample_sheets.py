import logging
from pathlib import Path
from typing import Union, List, Tuple

from sample_sheet import SampleSheet
from utils.seq_utils import get_rev_comp
from utils.utils import ExperimentalDesign

NUM_EXPECTED_SPACER = 4

MINIMUM_DIFF = 2

logger = logging.getLogger(__name__)


def test_for_collision(s1: str, s2: str, num_diff: int) -> bool:
    """
    Ensure sure that s1 and s2 differ by at more than <num_diff> nucleotides.
    :param s1: The first sequence.
    :param s2: The second sequence.
        len(s1) == len(s2)
    :param num_diff: Minimum number of nucleotide differences - 1.
    :return: Whether there are less than or equal to <num_diff> differences
        between <s1> and <s2>.
    """

    assert len(s1) == len(s2)
    num_eq = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            num_eq += 1
    return num_eq >= len(s1) - num_diff


def find_collisions(sample_sheet_path: Path) -> None:
    """
    Reports all instances where there are index collisions.
    :param sample_sheet_path: The path to the sample sheet to be checked for error.
    :return: None
    """
    sample_sheet = SampleSheet(sample_sheet_path)

    samples = sorted(sample_sheet.samples,
                     key=lambda x: x.sample_id.split('_')[1])
    sample_names = [sample.sample_id.split('_')[1]
                    for sample in samples]

    def get_collisions(indexes: List[str]) -> List[Tuple[int, int]]:
        collisions = []
        for i, index1 in enumerate(indexes):
            for j, index2 in enumerate(indexes[i + 1:]):
                if test_for_collision(index1, index2, MINIMUM_DIFF):
                    collisions.append((i, j))
        return collisions

    def format_collisions(indexes: List[str], samples: List[str],
                          collisions: List[Tuple[int, int]]) -> str:
        rtrn = ''
        for i, j in collisions:
            index1 = indexes[i]
            index2 = indexes[j]

            sample1 = samples[i]
            sample2 = samples[j]

            bars = ''
            for k, bp in enumerate(index1):
                if bp == index2[k]:
                    bars += '|'
                else:
                    bars += ' '

            s = ' '.join([sample1, '-', sample2,
                          '\n\t', index1,
                          '\n\t', bars,
                          '\n\t', index2, '\n\n'])

            rtrn += s

        return rtrn

    if sample_sheet.is_paired_end:
        i5_indices = [sample.index2 for sample in sample_sheet.samples]
        if get_collisions(i5_indices):
            print('=== i5 Collisions ===')
            print(format_collisions(i5_indices, sample_names,
                                    get_collisions(i5_indices)))

    i7_indices = [sample.index2 for sample in sample_sheet.samples]
    if get_collisions(i7_indices):
        print('=== i7 Collisions ===')
        print(format_collisions(i7_indices, sample_names,
                                get_collisions(i7_indices)))


def apply_basic_formatting(sample_sheet_path: Path, directory: Path = None,
                           verbose: bool = False) -> Path:
    """
    Performs several formatting steps to remove superfluous lines from the
    given sample sheet to allow for interpretation by the demultiplexer and
    sample_sheet module.

    1. Remove pre-header line.
        Excel often inserts a header line similar to the name of sample sheet.
        as the first line in the csv. This needs to be removed.

    2. Remove alternating spacer lines.
        Excel will often insert unnecessary a spacer line every other line.

    :param sample_sheet_path: Path to original sample sheet.
    :param directory: Path to the output directory. Parent of
        <sample_sheet_path> by default.
    :param verbose: Whether to print deatil about pruning process to console.
    :return: Path to the pruned sample sheet.
    """

    with open(sample_sheet_path, 'r') as sample_sheet:
        lines = sample_sheet.readlines()
    for i, line in enumerate(lines):
        lines[i] = line.strip()

    # === Remove header line ===
    while lines[0][:8] != "[Header]":
        msg = f"Pre-header line removed: {lines[0]}"
        if verbose:
            print(msg)
        logger.info(msg)
        lines.pop(0)

    # === Remove Purposeless Comma Rows ===
    def is_commas(s: str) -> bool:
        return all([c == "," or not c.strip()
                    for c in s.strip()])

    num_spacer_rows = 0
    for row in lines:
        if is_commas(row):
            num_spacer_rows += 1

    # If there are the correct number of spacer rows do not create new file.
    if num_spacer_rows <= NUM_EXPECTED_SPACER:
        msg = f"Found {num_spacer_rows} spacer rows. {sample_sheet_path.name}" \
              f" is properly formatted. No modification will be performed."
        if verbose:
            print(msg)
        logger.info(msg)
        return sample_sheet_path

    # First line that is all commas.
    first_spacer_ind = [is_commas(row) for row in lines].index(True)
    i = first_spacer_ind
    while is_commas(lines[i]):
        lines.pop(i)
        i += 1
    i -= first_spacer_ind

    while not lines[-1].strip():
        lines.pop()

    while is_commas(lines[-1]):
        lines.pop()

    # === Remove Purposeless Comma Columns ===

    def additional_cols() -> bool:
        # Are there tailing columns that are all commas
        for line in lines:
            if line[-1] != ',':
                return False
        return True

    while additional_cols():
        for i, line in enumerate(lines):
            lines[i] = line[:-1].strip()

    # === Create new csv ==

    if not directory:
        directory = sample_sheet_path.parent
    new_path = directory / (sample_sheet_path.name[:-4] +
                            '_trimmed.csv')

    msg = f"Sample sheet contained extra lines, removed {i} extra lines and " \
          f"created new sample sheet at {new_path.name} in the same directory."
    if verbose:
        print(msg)
    logger.info(msg)

    lines = [line + '\n' for line in lines]

    with open(new_path, 'w') as out_file:
        out_file.writelines(lines)

    return new_path


def reverse_compliment(sample_sheet_path: Path, directory: Union[str, Path],
                       rev_i5: bool, rev_i7: bool, verbose: bool = False) \
        -> Path:
    """
    Converts all i5 primers in the <sample_sheet> to their reverse compliment
    and saves the new sample sheet.

    :param rev_i7: Whether to reverse the i7 index.
    :param rev_i5: Whether to reverse the i5 index.
    :param sample_sheet_path: Path to sample sheet. Must be paired end.
    :param directory: Path to the output directory. Parent of
        <sample_sheet_path> by default.
    :param verbose: Whether to print info to console.
    :return: Path to newly create sample sheet.
    """

    sample_sheet = SampleSheet(sample_sheet_path)

    if not sample_sheet.is_paired_end:
        raise ValueError("Received a single end sample sheet.")

    if rev_i5:
        for sample in sample_sheet.samples:
            index_2 = sample.index2
            rc_index_2 = get_rev_comp(index_2)
            sample.index2 = rc_index_2

            msg = f"{sample.Sample_ID}\t{index_2}\t->\t{rc_index_2}"
            if verbose:
                print(msg)
            logger.debug(msg)

    if rev_i7:
        for sample in sample_sheet.samples:
            index_1 = sample.index1
            rc_index_1 = get_rev_comp(index_1)
            sample.index1 = rc_index_1

            msg = f"{sample.Sample_ID}\t{index_1}\t->\t{rc_index_1}"
            if verbose:
                print(msg)
            logger.debug(msg)

    if not directory:
        directory = sample_sheet_path.parent
    new_path = directory / (sample_sheet_path.name[:-4] +
                            '_rev_comp.csv')

    sample_sheet.write(open(new_path, 'w'))

    return new_path


def fix_sample_sheet(sample_sheet_path: Path,
                     rev_i5: bool = False, rev_i7: bool = False,
                     directory: Path = None,
                     verbose: bool = False) -> Path:
    """
    Performs several changes to the given sample sheet to prepare it for 
    use in demultiplexing.
    
    1. Remove pre-header line.
        Excel often inserts a header line similar to the name of sample sheet.
        as the first line in the csv. This needs to be removed.

    2. Remove alternating spacer lines.
        Excel will often insert unnecessary a spacer line every other line.

    3. Remove empty rows at the end of the file.

    4. Remove empty columns to the right of the data.
    
    :param rev_i7: Whether to convert i7 indices to their reverse compliments.
    :param rev_i5: Whether to convert i5 indices to their reverse compliments.
    :param sample_sheet_path: Path to the sample sheet to be modified.
    :param directory: Path to the output directory. Parent of
        <sample_sheet_path> by default.
    :param rev_comp: Whether to apply reverse compliment to index2 primers iif
        sample sheet is paired end.
    :param verbose: Whether to display progress.
    :return: Path to modified sample sheet.
    """
    pruned_sample_sheet = apply_basic_formatting(sample_sheet_path, directory,
                                                 verbose)

    if rev_i5 and not SampleSheet(pruned_sample_sheet).is_paired_end:
        raise ValueError("Cannot reverse i5 indexes on single ended "
                         "sample sheet.")
    if rev_i5 or rev_i7:
        pruned_sample_sheet = \
            reverse_compliment(pruned_sample_sheet,
                               rev_i5=rev_i5, rev_i7=rev_i7,
                               directory=directory,
                               verbose=verbose)

    find_collisions(pruned_sample_sheet)

    return pruned_sample_sheet

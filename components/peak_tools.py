from pathlib import Path
from typing import Tuple, List

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

from utils.job_manager import JobManager
from utils.path_manager import cmdify


def get_common_peaks(peaksfiles: List[Path], out_peaksfile: Path) \
        -> Tuple[Path, str]:

    cmd = cmdify(
        'cat', *peaksfiles, '|',
        'bedtools sort '
    )

    for peaksfile in peaksfiles:
        cmd += cmdify(
            ' | bedtools intersect',
            '-a stdin',
            '-b', peaksfile, ' -wa',
            '| bedtools sort '
        )
    cmd += cmdify('| cut -f1-3 | bedtools sort | bedtools merge >', out_peaksfile)
    return out_peaksfile, cmd


def get_merged_peaks(peaksfiles: List[Path], out_peaksfile: Path):
    cmd = cmdify(
        'cat', *peaksfiles, '| cut -f1-3 | bedtools sort | bedtools merge >', out_peaksfile
    )
    return out_peaksfile, cmd


def rem_ext(filepath: Path) -> str:
    return filepath.name.split('.')[0]


def get_venn_peaks(peakfile_a: Path, peakfile_b: Path, title: str,
                   out_dir: Path, jobs_manager: JobManager,
                   set_a_name: str = None, set_b_name: str = None) \
        -> Tuple[Path, Path, Path]:
    """
    Generates a Venn diagram comparing the peaks in two files and outputs files
    representing unique and common peaks.

    This function compares two peak files (typically in BED format) using bedtools,
    identifies unique peaks in each file, and the peaks common to both. It then
    generates a Venn diagram representing these relationships and saves the diagram
    as a PNG file. Additionally, it outputs the unique and common peak files.

    Parameters:
    :param set_a_name:
    :param peakfile_a: (Path) Path to the first peak file.
    :param peakfile_b: (Path) Path to the second peak file.
    :param title: (str) Title for the Venn diagram and part of the filename for
     output files.
    :param out_dir: (Path) Directory path where the output files will be saved.
    :param set_a_name: (str, optional) Name for the first set in the Venn diagram.
                                Defaults to the name of peakfile_a without its extension.
    :param set_b_name: (str, optional) Name for the second set in the Venn diagram.
                                Defaults to the name of peakfile_b without its extension.

    :returns Tuple[Path, Path, Path]: Paths to the output files for peaks
    unique to file A, peaks unique to file B, and peaks common to both
    files, respectively.

    Note:
    The function requires bedtools to be installed and accessible from the command line.
    It also relies on matplotlib for generating the Venn diagram.

    Example usage:
    >>> a_only, b_only, common = get_venn_peaks(Path("fileA.bed"), Path("fileB.bed"), "Comparison of A and B", Path("/output/directory"))
    """

    cmds = []

    if set_a_name is None:
        set_a_name = rem_ext(peakfile_a)

    if set_b_name is None:
        set_b_name = rem_ext(peakfile_b)

    out_dir = out_dir / title.replace(' ', '_')
    cmds.append(cmdify('mkdir', out_dir))

    a_only = out_dir / ('only_' + peakfile_a.name)
    b_only = out_dir / ('only_' + peakfile_b.name)
    common = out_dir / f"common_{rem_ext(peakfile_a)}_{rem_ext(peakfile_b)}.bed"

    cmds.append(cmdify('bedtools intersect',
            '-a', peakfile_a,
            '-b', peakfile_b,
            '-wa', '-v',
            '| cut -f1-3 >', a_only))

    cmds.append(cmdify('bedtools intersect',
            '-a', peakfile_b,
            '-b', peakfile_a,
            '-wa', '-v',
            '| cut -f1-3 >', b_only))

    cmds.append(cmdify(
        'cat', peakfile_a, peakfile_b,
        '| bedtools sort',
        '| bedtools intersect',
        '-a stdin',
        '-b', peakfile_a,
        '-wa',
        '| bedtools intersect',
        '-a stdin',
        '-b', peakfile_b,
        '-wa',
        '| bedtools sort',
        '| bedtools merge',
        '| cut -f1-3 >', common))

    for cmd in cmds:
        jobs_manager.execute_lazy(cmd)

    def get_num_lines(filepath: Path):
        with open(filepath, 'r') as file:
            return len(file.readlines())
    num_a = get_num_lines(a_only)
    num_b = get_num_lines(b_only)
    num_common = get_num_lines(common)

    plt.figure(figsize=(8, 6))
    venn = venn2(subsets=(num_a, num_b, num_common),
                 set_labels=(set_a_name, set_b_name))
    plt.title(title)

    # Customizing the Venn diagram's appearance

    patch_a = venn.get_patch_by_id('10')
    if patch_a is not None:
        patch_a.set_color('skyblue')
        patch_a.set_edgecolor('black')

    # Attempt to retrieve and color the patch for set B
    patch_b = venn.get_patch_by_id('01')
    if patch_b is not None:
        patch_b.set_color('pink')
        patch_b.set_edgecolor('black')

    # Attempt to retrieve and color the patch for the intersection
    intersection_patch = venn.get_patch_by_id('11')
    if intersection_patch is not None:
        intersection_patch.set_color('plum')
        intersection_patch.set_alpha(0.7)
        # Set transparency for intersection

    # Adjust the font size for labels and title
    if venn.set_labels:
        for text in filter(None,
                           venn.set_labels):  # This filters out None values
            text.set_fontsize(14)

    # Adjust the font size for subset labels
    if venn.subset_labels:
        for text in filter(None,
                           venn.subset_labels):  # This filters out None values
            text.set_fontsize(12)

    plt.title(title, fontsize=16)

    # Save and show the improved Venn diagram
    plt.savefig(out_dir / (title.replace(' ', '_') + '.png'),
                bbox_inches='tight')  # Save with a tight layout
    plt.show()
    plt.clf()

    return a_only, common, b_only



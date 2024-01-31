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
    cmd += cmdify(' | bedtools sort | bedtools merge >', out_peaksfile)
    return out_peaksfile, cmd

def get_merged_peaks(peaksfiles: List[Path], out_peaksfile: Path):
    cmd = cmdify(
        'cat', *peaksfiles, '| bedtools sort | bedtools merge >', out_peaksfile
    )
    return out_peaksfile, cmd

def rem_ext(filepath: Path) -> str:
    return filepath.name.split('.')[0]

def get_venn_peaks(peakfile_a: Path, peakfile_b: Path, title: str,
                   out_dir: Path, jobs_manager: JobManager,
                   set_a_name: str = None, set_b_name: str = None) \
        -> Tuple[Path, Path, Path, List[str]]:
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
    >>> a_only, b_only, common, cmds = get_venn_peaks(Path("fileA.bed"), Path("fileB.bed"), "Comparison of A and B", Path("/output/directory"))
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

    cmds += cmdify('bedtools intersect',
            '-a', peakfile_a,
            '-b', peakfile_b,
            '-wa', '-v',
            '>', a_only)

    cmds += cmdify('bedtools intersect',
            '-a', peakfile_b,
            '-b', peakfile_a,
            '-wa', '-v',
            '>', b_only)

    cmds += cmdify(
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
        '>', common)

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
    venn.get_patch_by_id('10').set_color('skyblue')  # Set color for set A
    venn.get_patch_by_id('10').set_edgecolor('black')
    venn.get_patch_by_id('01').set_color('pink')  # Set color for set B
    venn.get_patch_by_id('01').set_edgecolor('black')
    venn.get_patch_by_id('11').set_color('plum')  # Set color for intersection
    venn.get_patch_by_id('11').set_alpha(
        0.7)  # Set transparency for intersection

    # Adjust the font size for labels and title
    for text in venn.set_labels:
        text.set_fontsize(14)
    for text in venn.subset_labels:
        text.set_fontsize(12)

    plt.title(title, fontsize=16)

    # Save and show the improved Venn diagram
    plt.savefig(out_dir / (title.replace(' ', '_') + '.png'),
                bbox_inches='tight')  # Save with a tight layout
    plt.show()
    plt.clf()

    return a_only, common, b_only, cmds



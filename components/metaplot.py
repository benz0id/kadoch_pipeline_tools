import logging
import sys
from pathlib import Path
from threading import Thread
from time import sleep
from typing import Dict, Callable, Union, Tuple, List

from components.deeptools import generate_bed_matrix
from utils.fetch_files import get_matching_files
from utils.job_formatter import JobBuilder, PythonJob, ExecParams
from utils.job_manager import JobManager
from utils.mark_colour_map import get_colour
from utils.path_manager import cmdify, PathManager
from utils.utils import TargetedDesign

logger = logging.getLogger(__file__)


def add(dic, key, val):
    if key in dic:
        dic[key].append(val)
    else:
        dic[key] = [val]


def heatmaps_by_mark(design: TargetedDesign,
                     mark_to_bigwigs_dir: Dict[str, Path],
                     mark_to_peaks: Dict[str, Union[Path,
                                              List[Tuple[str, Path]]]],
                     mark_to_scale_factor: Dict[str, int],
                     out_dir: Path, job_manager: JobManager,
                     job_builder: JobBuilder, file_manager: PathManager,
                     start_array: Callable,
                     stop_array: Callable,
                     cores_per_job: int = 1,
                     verbose: bool = False,
                     namesafe_check: bool = True) -> None:
    """
    Generate heatmaps for each mark.

    :param cores_per_job:
    :param design: The experimental design to be used.
    :param mark_to_bigwigs_dir: Maps a mark to the directory containing the
        bigwigs to use for that mark.
    :param mark_to_peaks:
        Maps a mark to a directory containing bed files or a list of
        bedfiles alongside their row names.

        Iff the dict maps to a path, extract bedfiles for each treatment for
        the given mark and merge them into a single common peaks bedfile.

        Iff the dict maps to a list of Tuple[str, Path], use each of the Paths
        as the row bedfile and the name as the name of that column.

    :param mark_to_scale_factor: Provides a scaling factor for each of the
        given marks.
    :param out_dir: The directory into which all the output figures and
        intermediary files will be placed.
    :param job_manager: The job manager to use when executing jobs.
    :param job_builder: The job manager to use when building jobs.
    :param file_manager: Used to manage and navigate the project directory.
    :param start_array: Command to call in order to start a job array.
    :param stop_array: Command to call in oder to stop a job array.
    :param verbose: Print addiitonal information regarding the process to the
        console.
    :param namesafe_check: Perform addition checks to ensure that the correct
        files are being used independent of specifications provided by
        <design>.
    :return: None.
    """
    namesafe_violations = {}

    def get_merged_peaks(directory, mark) -> Tuple[Path, str, str]:
        # Merge all the peaks files associated with <mark> into a single
        # peakfile.
        sample_names = design.query(get='sample_name',
                                    filters={'target': mark})

        treatments = [design.query(get='treatment',
                                   filters={'sample_name': sn})[0]
                      for sn in sample_names]

        peaks = get_matching_files(directory, sample_names,
                                   filetype='(((narrow|broad)Peak)|peaks.bed)',
                                   under_delim=True, paths=True,
                                   one_to_one=True)

        merged_peaks = out_dir / mark / 'common_peaks.bed'

        job_manager.execute_lazy(
            cmdify('cat', *peaks,
                   '| bedtools sort',
                   '| bedtools merge',
                   '>', merged_peaks))

        with open(merged_peaks) as peaksfile:
            n_merged_peaks = len(peaksfile.readlines())
        row_name = f'All Merged {mark} Peaks (n = {n_merged_peaks})'

        for i, bedfile in enumerate(peaks):
            components = bedfile.name.split('_')
            treatment = treatments[i]
            if treatment not in components or mark not in components:
                add(namesafe_violations, mark, (bedfile, treatment))

        # Display merging process.
        merge_str = ''
        merge_str += '\t --- Peak Merging ---\n'
        for file in peaks:
            with open(file, 'r') as peakfile:
                sub_n_peaks = len(peakfile.readlines())
            merge_str += '\t' + str(sub_n_peaks) + '\t' + file.name + '\n'
        merge_str += f'Final Merged:\n\t{n_merged_peaks}\t' \
                     f'{merged_peaks.name}\n'
        if verbose:
            print(merge_str)
        logger.info(merge_str)

        row_desc = f'Rows:\n\n' \
            f'\t{row_name}\t-\t{merged_peaks.name}\n\n'

        # Return path to merged peaks, name of merged peaks, and a description
        # to be used as a row name.
        return merged_peaks, row_name, row_desc

    def get_annotated_peaks(peaks: List[Tuple[str, Path]]) -> \
        Tuple[List[Path], List[str], str]:

        bedfiles = []
        row_names = []
        row_desc = f'Rows:\n\n'

        for rowname, bedfile in peaks:
            bedfiles.append(bedfile)

            with open(bedfile) as peaksfile:
                n_merged_peaks = len(peaksfile.readlines())
            row_name = f'{rowname} (n = {n_merged_peaks})'
            row_names.append(row_name)
            row_desc += f'\t{row_name}\t-\t{bedfile.name}\n\n'

        return bedfiles, row_names, row_desc


    for mark in set(design.get_marks()):
        assert mark in mark_to_bigwigs_dir
        assert mark in mark_to_peaks

    threads = []

    cmds = []

    # Generate matrices.
    for mark in set(design.get_marks()):
        # Configure output directory
        mark_common_dir = file_manager.make(out_dir /
                                            f'common_peaks_sub_{mark}')

        # Fetch Peaks and bigwigs.
        bigwig_dir = mark_to_bigwigs_dir[mark]
        peaks_dir = mark_to_peaks[mark]

        sample_names = design.query(get='sample_name',
                                    filters={'target': mark})
        bigwigs = get_matching_files(bigwig_dir, sample_names, filetype='bw',
                                     under_delim=True, paths=True,
                                     one_to_one=True)

        treatments = [design.query(get='treatment',
                                   filters={'sample_name': sn})[0]
                      for sn in sample_names]

        # Handle creation of Rows
        peaks = mark_to_peaks[mark]
        if isinstance(peaks, Path):
            peaks, row_names, row_desc = get_merged_peaks(peaks, mark)
            row_names = [row_names]
            peaks = [peaks]
        elif isinstance(peaks, list):
            peaks, row_names, row_desc = get_annotated_peaks(peaks)
        else:
            raise ValueError('mark_to_peaks contains invalid dictionary '
                             'values.')

        # Make matrix.
        matrix = mark_common_dir / 'compute_matrix.gz'
        pj = PythonJob('make ' + str(matrix),
                       [],
                       generate_bed_matrix,

                       beds=peaks,
                       bigwigs=bigwigs,
                       column_names=treatments,
                       row_names=row_names,
                       out_path=matrix,
                       jobs=job_manager,
                       builder=job_builder,
                       start_array=start_array,
                       stop_array=stop_array,
                       path_manager=file_manager,
                       num_cores=cores_per_job
                       )

        threads.append(Thread(target=job_manager.execute_lazy, args=[pj]))

        # Make string representation of the figure.
        col_strs = ''
        for i, bw in enumerate(bigwigs):
            filename = bw.name
            col_name = treatments[i]
            col_strs += '\n\n\t' + col_name + '\t' + '-' + '\t' + filename

        s = f'\t===   Figure Configuration for {mark}   ===\n' \
            f'\n' \
            f'Columns:' \
            f'{col_strs}\n' \
            + row_desc

        if verbose:
            print(s)
        logger.info(s)

        # Perform check to ensure that correct labels and files have been used.
        for i, bigwig in enumerate(bigwigs):
            components = bigwig.name.split('_')
            treatment = treatments[i]
            if treatment not in components or mark not in components:
                add(namesafe_violations, mark, (bigwig, treatment))

    if namesafe_check and namesafe_violations:
        violation_string = '\t===  Namesafe Check Failed ===\n'
        for mark in namesafe_violations:
            violation_string += f'Failures for {mark}\n'
            for bw, treatment in namesafe_violations[mark]:
                violation_string += f'\t{treatment} or {mark} ' \
                                    f'not in {bw.name}\n'
        if verbose:
            print(violation_string, file=sys.stderr)
        logger.info(violation_string)

    sleep(1)

    for cmd in cmds:
        job_manager.execute_lazy(cmd)

    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    arrayable = \
        ExecParams(max_runtime=(0, 1, 0), num_cores=cores_per_job,
                   ram_per_core=1024 * 4, builder=job_builder, wait=False)

    # Now, create figures.
    for mark in set(design.get_marks()):
        mark_common_dir = file_manager.make(out_dir /
                                            f'common_peaks_sub_{mark}')
        matrix = mark_common_dir / 'compute_matrix.gz'
        img = mark_common_dir / 'metaplot.png'

        scale_factor = mark_to_scale_factor[mark]

        colour = 'white,' + get_colour(mark)

        cmd = cmdify(
            "plotHeatmap",
            "-m", matrix,
            "-o", img,
            "--sortRegions descend",
            "--sortUsing mean",
            "--yMin 0 --zMin 0",
            "--yMax", scale_factor,
            "--zMax", scale_factor,
            "--colorList", f'"{colour}"'
        )
        job_manager.execute_lazy(cmd, arrayable)





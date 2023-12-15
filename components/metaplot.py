import logging
import os
from pathlib import Path
from threading import Thread
from typing import Dict, Callable

from components.deeptools import generate_bed_matrix
from utils.fetch_files import get_matching_files
from utils.job_formatter import JobBuilder, PythonJob, ExecParams
from utils.job_manager import JobManager
from utils.mark_colour_map import get_colour
from utils.path_manager import cmdify, PathManager
from utils.utils import ExperimentalDesign, TargetedDesign

logger = logging.getLogger(__file__)


def heatmaps_by_mark(design: TargetedDesign,
                     mark_to_bigwigs_dir: Dict[str, Path],
                     mark_to_peaks: Dict[str, Path],
                     mark_to_scale_factor: Dict[str, int],
                     out_dir: Path, job_manager: JobManager,
                     job_builder: JobBuilder, file_manager: PathManager,
                     start_array: Callable,
                     stop_array: Callable,
                     cores_per_job: int = 1,
                     verbose: bool = False) -> None:
    """
    Generate heatmaps for each mark.
    :param cores_per_job:
    :param design: The experimental design to be used.
    :param mark_to_bigwigs_dir: Maps a mark to the directory containing the
        bigwigs to use for that mark.
    :param mark_to_peaks: Maps a mark to the directory containing the
        peaks to use for that mark (assumed to be either bed, narrowPeak or
        broadPeak files)
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
    :return: None.
    """

    for mark in set(design.get_marks()):
        assert mark in mark_to_bigwigs_dir
        assert mark in mark_to_peaks

    threads = []

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
        peaks = get_matching_files(peaks_dir, sample_names,
                                   filetype='(((narrow|broad)Peak)|peaks.bed)',
                                   under_delim=True, paths=True,
                                   one_to_one=True)

        # Merge peaks.
        merged_peaks = mark_common_dir / 'common_peaks.bed'
        job_manager.execute_lazy(cmdify('cat', *peaks, '>', merged_peaks))
        with open(merged_peaks) as peaksfile:
            n_peaks = len(peaksfile.readlines())

        # Configure row and column names.
        row_name = f'All Merged {mark} Peaks (n = {n_peaks})'
        treatments = [design.query(get='treatment',
                                   filters={'sample_name': sn})[0]
                      for sn in sample_names]

        # Make matrix.
        matrix = mark_common_dir / 'compute_matrix.gz'
        pj = PythonJob('make ' + str(matrix),
                       [],
                       generate_bed_matrix,

                       beds=[merged_peaks],
                       bigwigs=bigwigs,
                       column_names=treatments,
                       row_names=[row_name],
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
            f'Rows:\n\n' \
            f'\t{row_name}\t-\t{merged_peaks}'

        if verbose:
            print(s)
        logger.info(s)

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





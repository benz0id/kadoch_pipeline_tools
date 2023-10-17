from pathlib import Path
from typing import List, Callable

from utils.fetch_files import get_unique_filename
from utils.job_formatter import JobBuilder, ExecParams, Runtime
from utils.job_manager import JobManager
from utils.path_manager import cmdify, PathManager


def generate_bed_matrix(beds: List[Path], bigwigs: List[Path],
                        column_names: List[str], row_names: List[str],
                        out_path: Path, jobs: JobManager, builder: JobBuilder,
                        start_array: Callable, stop_array: Callable,
                        path_manager: PathManager) -> None:
    """
    Quantify each of the given bigwig files against the given beds, to form a
    large matrix of counts for use in heatmap plotting.

    if beds = [bed1, bed2, bed3] and bigwigs = [wig1, wig2, wig3], then the
    final matrix will be of the form.

            wig1    wig2    wig3
    bed1
    bed2           <vals>
    bed3

    Where each bed is expanded to all the sites it contains along the y
    axis.

    This uses Alex's code as a skeleton.

    :param row_names:
    :param beds:
    :param bigwigs:
    :param column_names:
    :param out_path:
    :param jobs:
    :param builder:
    :param start_array:
    :param stop_array:
    :return:
    """

    four_core = ExecParams(max_runtime=Runtime(0, 0, 30),
                           num_cores=4,
                           ram_per_core=512,
                           builder=builder)

    to_remove = []

    matrix_files = []

    # Calculate submatrices.
    start_array()
    for i, bedfile in enumerate(beds):
        matrix_files.append([])
        for j, bigwig in enumerate(bigwigs):
            tmp = path_manager.purgeable_files_dir / (bedfile.name[:-4] + '_'
                                                      + bigwig.name[:-3] + '.gz')
            to_remove.append(tmp)

            cmd = cmdify(
                "computeMatrix",
                "reference-point",
                "--referencePoint", "center",
                "-b", "1000",
                "-a", "1000",
                "-p", four_core.num_cores,
                "-R", bedfile,
                "-S", bigwig,
                "--binSize 50",
                "--sortRegions", "keep",
                "--missingDataAsZero",
                "-o", tmp)
            cmd += cmdify(
                "computeMatrixOperations relabel"
                '-m', tmp,
                '--groupLabel', "'" + row_names[i] + "'",
                '--sampleLabel' "'" + column_names[j] + "'"
            )
            jobs.execute(cmd, four_core)
            matrix_files[i].append(tmp)
    stop_array()

    # Combine columns and rename.
    cols = []
    for i in range(len(matrix_files[0])):
        col = [row[i] for row in matrix_files]

        tmp_col = path_manager.purgeable_files_dir / (
                    get_unique_filename() + '.gz')

        cmd = cmdify(
            "computeMatrixOperations rbind",
            '-m', *col,
            '-o', tmp_col
        )
        jobs.execute(cmd)
        cols.append(tmp_col)

        to_remove.append(tmp_col)

    # Combine rows.
    cmd = cmdify(
        "computeMatrixOperations cbind",
        '-m', *cols,
        '-o', out_path
    )
    jobs.execute(cmd)


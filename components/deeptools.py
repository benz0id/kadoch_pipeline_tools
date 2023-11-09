from pathlib import Path
from typing import List, Callable

from utils.fetch_files import get_unique_filename, outpath_to_dirname
from utils.job_formatter import JobBuilder, ExecParams, Runtime
from utils.job_manager import JobManager
from utils.path_manager import cmdify, PathManager


def generate_bed_matrix(beds: List[Path], bigwigs: List[Path],
                        column_names: List[str], row_names: List[str],
                        out_path: Path, jobs: JobManager, builder: JobBuilder,
                        start_array: Callable, stop_array: Callable,
                        path_manager: PathManager, num_cores: int = 4) -> None:
    """
    Quantify each of the given bigwig files against the given beds, to form a
    large matrix of counts for use in heatmap plotting.

    if beds = [bed1, bed2, bed3] and bigwigs = [wig1, wig2, wig3], then the
    final matrix will be of the form.

            wig1    wig2    wig3
    bed1
    bed2           <vals>
    bed3

    Where each bed is expanded to all the sites it contains along the yaxis.

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
                           num_cores=num_cores,
                           ram_per_core=512,
                           builder=builder)

    one_core = ExecParams(max_runtime=Runtime(0, 0, 30),
                          num_cores=1,
                          ram_per_core=1048,
                          builder=builder)

    matrix_files = []

    tmp_dir = path_manager.purgeable_files_dir / outpath_to_dirname(out_path)

    # Calculate submatrices.
    start_array()
    for i, bedfile in enumerate(beds):
        matrix_files.append([])
        for j, bigwig in enumerate(bigwigs):
            tmp = tmp_dir / f'{j}.{i}.tmp.gz'
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
            cmd += '\n' + cmdify(
                "computeMatrixOperations relabel",
                '-m', tmp,
                '--groupLabel', "'" + row_names[i] + "'",
                '--sampleLabel', "'" + column_names[j] + "'",
                '-o', tmp
            )
            jobs.execute(cmd, four_core)
            matrix_files[i].append(tmp)
    stop_array()

    # Combine columns and rename.
    cols = []
    start_array()
    for i in range(len(matrix_files[0])):
        col = [row[i] for row in matrix_files]

        # There
        if len(col) == 1:
            cols.append(col[0])
            continue

        tmp_col = tmp_dir / (f'col_{i}.gz')

        cmd = cmdify(
            "computeMatrixOperations rbind",
            '-m', *col,
            '-o', tmp_col
        )
        jobs.execute(cmd, one_core)
        cols.append(tmp_col)

    stop_array()

    # Combine rows.
    start_array()
    cmd = cmdify(
        "computeMatrixOperations cbind",
        '-m', *cols,
        '-o', out_path
    )
    stop_array()
    jobs.execute(cmd, one_core)

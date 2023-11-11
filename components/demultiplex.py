from pathlib import Path
from typing import List, Callable

from utils.fetch_files import get_matching_files
from utils.job_formatter import ExecParams, JobBuilder
from utils.job_manager import JobManager
from utils.path_manager import PathManager
from wrappers.bcl2fastq import Demultiplexer
from wrappers.fastqc import FastQC
from utils.utils import ExperimentalDesign, get_model_from_sample_sheet, combine_runs 


def demult_and_fastqc(sample_sheet: Path, sequencing_results: Path,
                      builder: JobBuilder, path_manager: PathManager,
                      jobs_manager: JobManager, start_array: Callable,
                      drop_array: Callable) -> List[Path]:
    """
    Run generic demultiplexing and fastqc pipeline. Return paths to the
    resultant fastqs.

    === Project Structure Assumptions ===

    Fastqs in the fastqs directory are produce by this sequencing run.


    :param sample_sheet: Path to the sample sheet to use for demultiplexing.
    :param sequencing_results: Path to sequencing results directory.
    :param builder: The builder to use to build this component's jobs.
    :param path_manager: Path manager for accessing project directories.
    :param jobs_manager: Job manager to use when queueing jobs.
    :param start_array: Command to call to signify that a job array is
        starting.
    :param drop_array: Command to call to signify that the program can continue
        without completing the jobs currently in the array.
    :return: A list of fastqfiles produces by demultiplexing. Directory must 
        only contain raw output fastqs if this return is to be valid.
    """

    design = get_model_from_sample_sheet(sample_sheet)

    demult = Demultiplexer(path_manager)
    heavy_job = ExecParams(max_runtime=(0, 1, 0), num_cores=16,
                           ram_per_core=1024 * 4, builder=builder, wait=True)
    qc_job = ExecParams(max_runtime=(0, 1, 0), num_cores=4, ram_per_core=300,
                        builder=builder, wait=False)

    # === Demultiplex ===
    cmd = demult.get_demultiplex_cmd(
        sequencing_results, sample_sheet,
        path_manager.fastqs_dir, heavy_job.num_cores,
        reports_dir=path_manager.demult_stats / 'Reports',
        stats_dir=path_manager.demult_stats / 'Stats')
    jobs_manager.execute_lazy(cmd, heavy_job)

    # === FastQC ===
    fastqs = get_matching_files(path_manager.fastqs_dir, design.get_samples(),
                                ['Undetermined'], containing=True, paths=True)

    qc = FastQC()
    jobs = qc.fast_qc_arrayed(fastqs,
                              path_manager.fastqc_dir,
                              qc_job)

    start_array()
    for job in jobs:
        jobs_manager.execute_lazy(job, qc_job)
    drop_array()

    return fastqs

from pathlib import Path
from typing import List, Callable

from utils.fetch_files import get_matching_files
from utils.job_formatter import ExecParams
from utils.job_manager import JobManager
from utils.path_manager import PathManager
from wrappers.bcl2fastq import Demultiplexer
from wrappers.fastqc import FastQC


def demult_and_fastqc(sample_sheet: Path,
                      path_manager: PathManager,  jobs_manager: JobManager,
                      heavy_job: ExecParams, parallel_job: ExecParams,
                      start_array: Callable, drop_array: Callable) -> List[Path]:
    """
    Run generic demultiplexing and fastqc pipeline. Return paths to the
    resultant fastqs.

    === Project Structure Assumptions ===

    Fastqs in the fastqs directory are produce by this sequencing run.


    :param sample_sheet: Path to the sample sheet to use for demultiplexing.
    :param path_manager: Path manager for accessing project directories.
    :param jobs_manager: Job manager to use when queueing jobs.
    :param heavy_job: Execution parameters to use when executing heavy jobs.
        (demultiplexing)
    :param parallel_job: Execution paramers for mid-size job to be used in
        parallel.
    :param start_array: Command to call to signify that a job array is
        starting.
    :param drop_array: Command to call to signify that the program can continue
        without completing the jobs currently in the array.
    :return: A list of fastqfiles produces by demultiplexing.
    """
    demult = Demultiplexer(path_manager)

    # === Demultiplex ===
    cmd = demult.get_demultiplex_cmd(path_manager.sequencing_dir,
                                     sample_sheet,
                                     path_manager.fastqs_dir,
                                     heavy_job.num_cores)
    jobs_manager.execute_lazy(cmd, heavy_job)

    # === FastQC ===
    fastqs = get_matching_files(path_manager.fastqs_dir, ['fastq'],
                                ['Undetermined'], containing=True, paths=True)
    qc = FastQC()
    jobs = qc.fast_qc_arrayed(fastqs,
                              path_manager.fastqs_dir,
                              parallel_job)

    start_array()
    for job in jobs:
        jobs_manager.execute_lazy(job, parallel_job)
    drop_array()

    return fastqs

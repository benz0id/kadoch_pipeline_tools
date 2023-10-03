from pathlib import Path
from typing import List, Callable

from utils.job_formatter import ExecParams
from utils.job_manager import JobManager
from utils.path_manager import PathManager
from wrappers.bcl2fastq import Demultiplexer


def demult_and_fastqc(path_manager: PathManager,
                      jobs_manager: JobManager,
                      heavy_job: ExecParams, parallel_jobs: ExecParams,
                      start_array: Callable, drop_array: Callable) -> List[Path]:
    demult = Demultiplexer(path_manager)

    # === Demultiplex ===
    cmd = demult.get_demultiplex_cmd(path_manager.sequencing_dir,
                                     sample_sheet, path_manager.fastqs_dir,
                                     heavy_job.num_cores), \
        heavy_job
    lazy(demult.get_demultiplex_cmd(seq_res_dir, sample_sheet, fastqs_dir, arrayable_heavy_o2.num_cores), arrayable_heavy_o2)


    # === FastQC ===

    fastqc_dir = pm.sequencing_dir / 'fastqc'
    lazy(c('mkdir', fastqc_dir))

    fastqs = get_matching_files(fastqs_dir, ['fastq'], ['Undetermined'], containing=True, paths=True)

    qc_job = ExecParams(max_runtime=(0, 1, 0), num_cores=4, ram_per_core=300, builder=slurm, wait=False)
    qc = FastQC()


    jobs = qc.fast_qc_arrayed(fastqs,
                              fastqc_dir,
                              qc_job)

    slurm.begin_array()
    for job in jobs:
        lazy(job, qc_job)
    slurm.wait_for_all_jobs()

    exit(0)
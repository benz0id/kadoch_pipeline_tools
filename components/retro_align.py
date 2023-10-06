from pathlib import Path
from random import random
from typing import List, Callable

from utils.fetch_files import copy_to_cmds, get_matching_files
from utils.job_formatter import JobBuilder, ExecParams
from utils.job_manager import JobManager
from utils.path_manager import cmdify
from utils.utils import combine_cmds


def retro_align(fastqs: List[Path], alignment_dir: Path, cmd: str,
                aligments_results_dir: Path, jobs: JobManager,
                start_array: Callable, wait_array: Callable,
                job_builder: JobBuilder):

    bams_path = aligments_results_dir / 'bams'
    beds_path = aligments_results_dir / 'beds'
    bigwigs_path = aligments_results_dir / 'bigwigs'
    align_stats_path = aligments_results_dir / 'stats'

    jobs.execute_lazy(cmdify('mkdir', aligments_results_dir, bams_path, bigwigs_path, beds_path,
           align_stats_path))

    # Extract sample names from fastqs.
    sample_names = []
    for fastq in fastqs:
        name = fastq.name
        sample_name = name.split('_')[1]
        sample_names.append(sample_name)
    sample_names = sorted(sample_names)

    # Copy all matching sequencing result folders over.
    cmds = []
    cmds.extend(copy_to_cmds(bams_path,
                             get_matching_files(alignment_dir / 'bam',
                                                sample_names,
                                                containing=True, paths=True),
                             avoid_recopy=True))
    cmds.extend(copy_to_cmds(beds_path,
                             get_matching_files(alignment_dir / 'bed',
                                                sample_names,
                                                containing=True, paths=True),
                             avoid_recopy=True))
    cmds.extend(copy_to_cmds(bigwigs_path,
                             get_matching_files(alignment_dir / 'bw',
                                                sample_names,
                                                containing=True, paths=True),
                             avoid_recopy=True))
    cmds.extend(copy_to_cmds(align_stats_path,
                             get_matching_files(alignment_dir / 'stats_storage*',
                                                sample_names,
                                                containing=True, paths=True),
                             avoid_recopy=True))

    random.shuffle(cmds)
    combined_cmds = combine_cmds(cmds, 4)
    light_o2 = ExecParams(max_runtime=(0, 0, 10), num_cores=1,
                          ram_per_core=128, builder=job_builder)
    start_array()
    for cmd in combined_cmds:
        jobs.execute_lazy(cmd, light_o2)
    wait_array()

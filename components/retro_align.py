import collections
import os
from pathlib import Path
import random
from typing import List, Callable

from components.get_alignment_stats import get_alignment_stats
from utils.fetch_files import copy_to_cmds, get_matching_files
from utils.job_formatter import JobBuilder, ExecParams
from utils.job_manager import JobManager
from utils.path_manager import cmdify
from utils.utils import combine_cmds

AlignmentResults = collections.namedtuple('AlignmentResults', ['bam', 'bed',
                                                               'bw', 'stats'])


def retro_fetch_align_results(fastqs: List[Path], alignment_dir: Path,
                              aligments_results_dir: Path, jobs: JobManager,
                              start_array: Callable, wait_array: Callable,
                              job_builder: JobBuilder) -> AlignmentResults:

    bams_path = aligments_results_dir / 'bams'
    beds_path = aligments_results_dir / 'beds'
    bigwigs_path = aligments_results_dir / 'bigwigs'
    align_stats_path = aligments_results_dir / 'stats'

    res = AlignmentResults(bam=bams_path,
                           bed=beds_path,
                           bw=bigwigs_path,
                           stats=align_stats_path)

    jobs.execute_lazy(cmdify('mkdir', *res))

    # Extract sample names from fastqs.
    sample_ids = []
    for fastq in fastqs:
        name = fastq.name
        sample_name = '_'.join(name.split('_')[0:3])
        print(sample_name)
        sample_ids.append(sample_name)
    sample_ids = sorted(sample_ids)


    # Copy all matching sequencing result folders over.
    cmds = []
    cmds.extend(copy_to_cmds(bams_path,
                             get_matching_files(alignment_dir / 'bam',
                                                sample_ids,
                                                containing=True, paths=True),
                             avoid_recopy=True))
    cmds.extend(copy_to_cmds(beds_path,
                             get_matching_files(alignment_dir / 'beds',
                                                sample_ids,
                                                containing=True, paths=True),
                             avoid_recopy=True))
    cmds.extend(copy_to_cmds(bigwigs_path,
                             get_matching_files(alignment_dir / 'bw',
                                                sample_ids,
                                                containing=True, paths=True),
                             avoid_recopy=True))

    if 'atac' in alignment_dir.name:
        cmds.extend(copy_to_cmds(align_stats_path,
                                 get_matching_files(
                                     alignment_dir / 'stats_storage_atac',
                                     sample_ids,
                                     containing=True, paths=True),
                                 avoid_recopy=True))

    random.shuffle(cmds)
    combined_cmds = combine_cmds(cmds, 5)
    light_o2 = ExecParams(max_runtime=(0, 0, 10), num_cores=1,
                          ram_per_core=128, builder=job_builder)
    start_array()
    for cmd in combined_cmds:
        jobs.execute_lazy(cmd, light_o2)
    wait_array()

    get_alignment_stats([align_stats_path / align_stat for align_stat in
                         os.listdir(align_stats_path)],
                        out_filepath=aligments_results_dir / 'stats.txt')

    return res

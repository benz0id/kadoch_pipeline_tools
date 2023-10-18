import collections
import os
from pathlib import Path
import random
from typing import List, Callable, Tuple

from components.get_alignment_stats import get_alignment_stats
from utils.fetch_files import copy_to_cmds, get_matching_files
from utils.job_formatter import JobBuilder, ExecParams
from utils.job_manager import JobManager
from utils.path_manager import cmdify
from utils.utils import combine_cmds, ExperimentalDesign

AlignmentResults = collections.namedtuple('AlignmentResults',
                                          ['path', 'bam', 'bed', 'bw', 'stats'])


def make_alignment_results_directory(results_dir: Path) -> AlignmentResults:
    """
    Make a direcotry for storing alignment results.
    :param results_dir: Path to the alignment results dir.
    :return:
    """
    if not results_dir.exists():
        os.mkdir(results_dir)

    bams_path = results_dir / 'bam'
    beds_path = results_dir / 'bed'
    bigwigs_path = results_dir / 'bigwig'
    align_stats_path = results_dir / 'stat'

    res = AlignmentResults(path=results_dir,
                           bam=bams_path,
                           bed=beds_path,
                           bw=bigwigs_path,
                           stats=align_stats_path)

    for d in res:
        os.mkdir(d)


    return res


def retro_fetch_manual(design: ExperimentalDesign, out_dir: Path,
                       bam_dir: Path = None,
                       bai_dir: Path = None,
                       bed_dir: Path = None,
                       bigwig_dir: Path = None,
                       stats_dir: Path = None
                       ) -> Tuple[AlignmentResults, List[str]]:

    alignment_results = make_alignment_results_directory(out_dir)

    cmds = []

    def add_move_file_cmds(directory: Path, filetype: str, out_path: Path):
        end_regex = '.*' + filetype + '$'
        regexes = ['.*' + sample + end_regex
                   for sample in design.get_samples()]

        matching_files = get_matching_files(directory, regexes, paths=True)
        try:
            matching_files = design.align_to_samples(matching_files)
        except ValueError as e:
            print("Unable to find requested files.")
            get_matching_files(directory, regexes, paths=True, verbose=True)
            raise e

        for file in matching_files:
            cmds.append(cmdify('cp', file, out_path))

    if bam_dir:
        add_move_file_cmds(bam_dir, 'bam', alignment_results.bam)
    if bai_dir:
        add_move_file_cmds(bai_dir, 'bai', alignment_results.bam)
    if bed_dir:
        add_move_file_cmds(bed_dir, 'bed', alignment_results.bed)
    if bigwig_dir:
        add_move_file_cmds(bigwig_dir, 'bed', alignment_results.bw)
    if stats_dir:
        add_move_file_cmds(stats_dir, 'txt', alignment_results.stats)

    return alignment_results, cmds


def retro_fetch_align_results(fastqs: List[Path], alignment_dir: Path,
                              aligments_results_dir: Path, jobs: JobManager,
                              start_array: Callable, wait_array: Callable,
                              job_builder: JobBuilder) -> AlignmentResults:

    res = make_alignment_results_directory(results_dir=aligments_results_dir)

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
    cmds.extend(copy_to_cmds(res.bam,
                             get_matching_files(alignment_dir / 'bam',
                                                sample_ids,
                                                containing=True, paths=True),
                             avoid_recopy=True))
    cmds.extend(copy_to_cmds(res.bed,
                             get_matching_files(alignment_dir / 'beds',
                                                sample_ids,
                                                containing=True, paths=True),
                             avoid_recopy=True))
    cmds.extend(copy_to_cmds(res.bw,
                             get_matching_files(alignment_dir / 'bw',
                                                sample_ids,
                                                containing=True, paths=True),
                             avoid_recopy=True))

    if 'atac' in alignment_dir.name:
        cmds.extend(copy_to_cmds(res.stats,
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

    get_alignment_stats([res.stats / align_stat for align_stat in
                         os.listdir(res.stats)],
                        out_filepath=aligments_results_dir / 'stats.txt')

    return res

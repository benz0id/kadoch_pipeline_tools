from pathlib import Path

from utils.job_formatter import ExecParams, JobBuilder
from utils.fetch_files import get_matching_files, copy_to_cmds
from components.peak_calling import PeakCaller
from utils.path_manager import cmdify
from time import sleep
from components.get_alignment_stats import get_alignment_stats, get_sample_names
from wrappers.bcl2fastq import Demultiplexer
from wrappers.fastqc import FastQC
from components.merge_fastqs import merge_fastqs
from constants.old_storage_paths import *
from utils.utils import combine_cmds
from components.peaks_pca import PeakPCAAnalyser

import utils.cache_manager
import utils.path_manager
import utils.slurmify
import utils.job_manager
import os
import random

# Project Directory
wd = Path("/home/bet488/kadoch/ben/projects/grace_oncohis/atac")
raw_atac_dir = Path("/n/data1/dfci/pedonc/kadoch/atac_fastq")

# Pipeline Managers
pm = utils.path_manager.PathManager(wd, verbose=False)
cm = utils.cache_manager.CacheManager(pm, strict=False, name='atac_pipeline_cache') # TODO try strict
slurm = utils.slurmify.Slurmifier(pm, mailtype=["ALL"],
                               redirect_to_console=False)
jobs = utils.job_manager.JobManager(cm, pm, verbose=True)

# Job Configurations
generic_job = ExecParams(max_runtime=(0, 0, 10), num_cores=1, ram_per_core=128, builder=JobBuilder())
light_o2 = ExecParams(max_runtime=(0, 0, 10), num_cores=1, ram_per_core=128, builder=slurm)
arrayable_light_o2 = ExecParams(max_runtime=(0, 1, 0), num_cores=1, ram_per_core=128, builder=slurm, wait=False)
arrayable_heavy_o2 = ExecParams(max_runtime=(0, 1, 0), num_cores=8, ram_per_core=1024 * 4, builder=slurm, wait=False)

# Funciton aliases.
c = cmdify
run = jobs.execute
lazy = jobs.execute_lazy

# Experiment Meta-Info
conditions = ["empty", "D81WT", "D81A", "D81N"]

# Componenets
peak_caller = PeakCaller(jobs, pm, generic_job, arrayable_heavy_o2, lazy=True)
peak_counter = PeakPCAAnalyser(jobs, pm, generic_job, arrayable_light_o2, arrayable_heavy_o2,
                             start_job_array_fxn=slurm.begin_array,
                             wait_job_array_fxn=slurm.wait_for_all_jobs)

# === FastQC ===

fastqs_dir = wd / 'fastq'

fastqs = get_matching_files(fastqs_dir, ['fastq'], ['Undetermined'], containing=True, paths=True)

qc_job = ExecParams(max_runtime=(0, 1, 0), num_cores=4, ram_per_core=300, builder=slurm, wait=False)
qc = FastQC()

lazy(c('mkdir', 'fastqc'))

jobs = qc.fast_qc_arrayed(fastqs, 'fastqc', qc_job)

slurm.begin_array()
for job in jobs:
    lazy(job, qc_job)
slurm.wait_for_all_jobs()

# === Fetch Alignment Results ===

alignment_results_path = wd / 'alignment_results'
bams_path = alignment_results_path / 'bams'
beds_path = alignment_results_path / 'beds'
bigwigs_path = alignment_results_path / 'bigwigs'
align_stats_path = alignment_results_path / 'stats'

lazy(c('mkdir', alignment_results_path, bams_path, bigwigs_path, beds_path, align_stats_path))

# Extract sample names from fastqs.
sample_names = []
for fastq in fastqs:
    name = fastq
    sample_name = name.split('_')[1]
    sample_names.append(sample_name)
sample_names = sorted(sample_names)

# Copy all matching sequencing result folders over.
cmds = []
cmds.extend(copy_to_cmds(bams_path, get_matching_files(ATAC_BAMS, sample_names, containing=True, paths=True), avoid_recopy=True))
cmds.extend(copy_to_cmds(beds_path, get_matching_files(ATAC_BEDS, sample_names, containing=True, paths=True), avoid_recopy=True))
cmds.extend(copy_to_cmds(bigwigs_path, get_matching_files(ATAC_BIGWIGS, sample_names, containing=True, paths=True), avoid_recopy=True))
cmds.extend(copy_to_cmds(align_stats_path, get_matching_files(ATAC_STATS, sample_names, containing=True, paths=True), avoid_recopy=True))

random.shuffle(cmds)
combined_cmds = combine_cmds(cmds, 4)
slurm.begin_array()
for cmd in combined_cmds:
    lazy(cmd, arrayable_light_o2)
slurm.wait_for_all_jobs()

# === Begin Peak Calling ===

CUTOFF = 0.01
all_beds = []
all_bams = []

slurm.begin_array()
for cond in conditions:
    bams = get_matching_files(bams_path, [cond], ['bai', 'pdf', 'txt', 'Input'],
                              containing=True, paths=True, verbose=False)
    # print(' - '.join([cond, '\n\t' + '\n\t'.join(bams)]))
    beds = peak_caller.run_peak_calling(bams, cond + '_peaks',
                                        cutoff_q_val=CUTOFF, paired_end=True,
                                        peak_type='narrow', use_shift_model=False)
    all_beds.extend(beds)
    all_bams.extend(bams)
slurm.wait_for_all_jobs()

print('\n'.join([str(bed) for bed in all_beds]))

# === Get Counts Matrix ===

pca_dir = wd / 'pca_analysis'
peak_counter.do_peaks_pca_analysis(all_beds, all_bams, pca_dir)















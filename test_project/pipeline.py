from pathlib import Path

from components.merge_fastqs import merge_fastqs_grouped
from constants.data_paths import RC_DATA_CREDS
from utils.fetch_files import get_unique_filename
from utils.file_transfer import SSHFileTransferManager
from utils.job_formatter import ExecParams, JobBuilder

import utils.cache_manager
import utils.path_manager
import utils.slurmify
import utils.job_manager
from utils.utils import get_model_from_sample_sheet, combine_runs
from components.deeptools import generate_bed_matrix

# Project Directory
wd = Path("/Users/btudorpr/PycharmProjects/kadoch_pipeline_tools/test_project")

# Pipeline Managers
pm = utils.path_manager.PathManager(wd)
cm = utils.cache_manager.CacheManager(pm, strict=True) # TODO try strict
sl = utils.slurmify.Slurmifier(pm, mailtype=["ALL"],
                               redirect_to_console=False)
jobs = utils.job_manager.JobManager(cm, pm)

get_unique_filename()

# Job Configurations
generic_job = ExecParams(max_runtime=(0, 0, 10), num_cores=1, ram_per_core=128,
                         builder=JobBuilder())
light_o2 = ExecParams(max_runtime=(0, 0, 10), num_cores=1, ram_per_core=128,
                      builder=sl)
generic_o2 = ExecParams(max_runtime=(0, 1, 0), num_cores=1, ram_per_core=128,
                        builder=sl)

# Useful paths.
pd = pm.project_dir

# === Begin Pipeline ===

p = lambda *x: [Path(i) for i in x]

jobs.disable_execution()
generate_bed_matrix(
    beds=p('a.bed', 'b.bed', 'c.bed'),
    bigwigs=p('d.bw', 'e.bw', 'f.bw'),
    column_names=['g', 'h', 'i'],
    out_path=Path('a'),
    jobs=jobs,
    builder=JobBuilder(),
    start_array=lambda: print('s'),
    stop_array=lambda: print('w'),
    path_manager=pm
)




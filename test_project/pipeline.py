from pathlib import Path

from constants.data_paths import RC_DATA_CREDS
from utils.file_transfer import SSHFileTransferManager
from utils.job_formatter import ExecParams, JobBuilder

import utils.cache_manager
import utils.path_manager
import utils.slurmify
import utils.job_manager

# Project Directory
soft = Path("/home/bet488/soft")
wd = Path("/home/bet488/soft/kadoch_pipeline_tools/test_project")

# Pipeline Managers
pm = utils.path_manager.PathManager(wd)
cm = utils.cache_manager.CacheManager(pm, strict=True) # TODO try strict
sl = utils.slurmify.Slurmifier(pm, mailtype=["ALL"], redirect_to_console=False)
jobs = utils.job_manager.JobManager(cm, pm)

# Job Configurations
generic_job = ExecParams(max_runtime=(0, 0, 10), num_cores=1, ram_per_core=128, builder=JobBuilder())
generic_o2 = ExecParams(max_runtime=(0, 1, 0), num_cores=1, ram_per_core=128, builder=sl)
heavy_o2 = ExecParams(max_runtime=(0, 1, 0), num_cores=30, ram_per_core=8192, builder=sl)

# Useful paths.
pd = pm.project_dir

# === Begin Pipeline ===
# jobs.execute('ls', generic_job)

# jobs.execute('ls -a', generic_job)

# === Testing Data Transfer ===

"""
remote_files = SSHFileTransferManager(RC_DATA_CREDS)
ck_dir = Path('/labs/cklab')
basic_fromrd3 = remote_files.get_files_ignore_existing_job([ck_dir / 'transfer_to_o2.sh',
                                        ck_dir / 'peaks.list.txt'],
                                       pd / 'storage')

dir_to_r3 = remote_files.put_files_ignore_existing_job(pd / 'storage' / 'test_dir',
                                       ck_dir / 'ben')

r3_to_loc_dir = remote_files.get_files_ignore_existing_job(ck_dir / 'ben' / 'test_dir',
                                       pd / 'storage2')

jobs.execute_lazy(basic_fromrd3, generic_job)
jobs.execute_lazy(dir_to_r3, generic_job)
jobs.execute_lazy(r3_to_loc_dir, generic_job)

remote_files.close()
"""

# === Testing Job Caching ===

jobs.execute_lazy('echo a', generic_job)
jobs.execute_lazy('echo b', generic_job)

jobs.execute_lazy('echo New Step!', generic_job)

jobs.execute_lazy('echo c', generic_job)
jobs.execute_lazy('echo d', generic_job)

testfile = pd / 'test.txt'

jobs.execute_purgeable('touch ' + str(testfile), testfile, generic_job, lazy=False)

# We're doing slurm jobs now. Since my mac doesn't have slurm, let's disable their execution.

jobs.execute('echo schlurmski', generic_o2)
jobs.execute('echo schlurmski', heavy_o2)


sl.stop_readers()
cm.purge_data()






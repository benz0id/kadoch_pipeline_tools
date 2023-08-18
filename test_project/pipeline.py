from pathlib import Path
from utils.job_formatter import ExecParams, JobBuilder

import utils.cache_manager
import utils.path_manager
import utils.slurmify
import utils.job_manager

# Project Directory
wd = Path("/Users/btudorpr/PycharmProjects/kadoch_pipeline_tools/test_project")

# Pipeline Managers
pm = utils.path_manager.PathManager(wd)
cm = utils.cache_manager.CacheManager(pm, strict=True) # TODO try strict
sl = utils.slurmify.Slurmifier(pm, mailtype=["ALL"])
jobs = utils.job_manager.JobManager(cm, pm)

# Job Configurations
generic_job = ExecParams(max_runtime=(0, 0, 10), num_cores=1, ram_per_core=128, builder=JobBuilder())
generic_o2 = ExecParams(max_runtime=(0, 1, 0), num_cores=1, ram_per_core=128, builder=sl)

# Useful paths.
pd = pm.project_dir

# === Begin Pipeline ===
# jobs.execute('ls', generic_job)

# jobs.execute('ls -a', generic_job)

jobs.execute_lazy('echo a', generic_job)
jobs.execute_lazy('echo b', generic_job)

jobs.execute_lazy('echo New Step!', generic_job)

jobs.execute_lazy('echo c', generic_job)
jobs.execute_lazy('echo d', generic_job)

testfile = pd / 'test.txt'

jobs.execute_purgeable('touch ' + str(testfile), testfile, generic_job, lazy=False)

# We're doing slurm jobs now. Since my mac doesn't have slurm, let's disable their execution.

jobs.execute('echo schlurmski', generic_o2)





sl.stop_readers()
cm.purge_data()






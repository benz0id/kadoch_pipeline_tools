from pathlib import Path

from components.merge_fastqs import merge_fastqs_grouped
from constants.data_paths import RC_DATA_CREDS
from utils.file_transfer import SSHFileTransferManager
from utils.job_formatter import ExecParams, JobBuilder

import utils.cache_manager
import utils.path_manager
import utils.slurmify
import utils.job_manager
from utils.utils import get_model_from_sample_sheet, combine_runs

# Project Directory
wd = Path("/Users/btudorpr/PycharmProjects/kadoch_pipeline_tools/test_project")

# Pipeline Managers
pm = utils.path_manager.PathManager(wd)
cm = utils.cache_manager.CacheManager(pm, strict=True) # TODO try strict
sl = utils.slurmify.Slurmifier(pm, mailtype=["ALL"],
                               redirect_to_console=False)
jobs = utils.job_manager.JobManager(cm, pm)


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

fastqs = [
    "20230817_CGATAC178_H358_Rep1_S1_R1_001.fastq.gz",
    "20230817_CGATAC178_H358_Rep1_S1_R2_001.fastq.gz",
    "20230817_CGATAC179_H358_Rep3_S2_R1_001.fastq.gz",
    "20230817_CGATAC179_H358_Rep3_S2_R2_001.fastq.gz",
    "20230817_CGATAC180_H1373_Rep1_S3_R1_001.fastq.gz",
    "20230817_CGATAC180_H1373_Rep1_S3_R2_001.fastq.gz",
    "20230817_CGATAC181_H1373_Rep3_S4_R1_001.fastq.gz",
    "20230817_CGATAC181_H1373_Rep3_S4_R2_001.fastq.gz",
    "20230817_CGATAC182_HOP62_Rep2_S5_R1_001.fastq.gz",
    "20230817_CGATAC182_HOP62_Rep2_S5_R2_001.fastq.gz",
    "20230817_CGATAC183_HOP62_Rep3_S6_R1_001.fastq.gz",
    "20230817_CGATAC183_HOP62_Rep3_S6_R2_001.fastq.gz",
    "20230817_CGATAC184_H23_Rep1_S7_R1_001.fastq.gz",
    "20230817_CGATAC184_H23_Rep1_S7_R2_001.fastq.gz",
    "20230817_CGATAC185_H23_Rep3_S8_R1_001.fastq.gz",
    "20230817_CGATAC185_H23_Rep3_S8_R2_001.fastq.gz",
    "20230817_CGATAC186_Lu65_Rep1_S9_R1_001.fastq.gz",
    "20230817_CGATAC186_Lu65_Rep1_S9_R2_001.fastq.gz",
    "20230817_CGATAC187_Lu65_Rep3_S10_R1_001.fastq.gz",
    "20230817_CGATAC187_Lu65_Rep3_S10_R2_001.fastq.gz",
    "20230817_CGATAC188_H2030_Rep1_S11_R1_001.fastq.gz",
    "20230817_CGATAC188_H2030_Rep1_S11_R2_001.fastq.gz",
    "20230817_CGATAC189_H2030_Rep3_S12_R1_001.fastq.gz",
    "20230817_CGATAC189_H2030_Rep3_S12_R2_001.fastq.gz",
    "20230817_CGATAC190_H2122_Rep1_S13_R1_001.fastq.gz",
    "20230817_CGATAC190_H2122_Rep1_S13_R2_001.fastq.gz",
    "20230817_CGATAC191_H2122_Rep3_S14_R1_001.fastq.gz",
    "20230817_CGATAC191_H2122_Rep3_S14_R2_001.fastq.gz",
    "20230817_CGATAC192_HCC44_Rep1_S15_R1_001.fastq.gz",
    "20230817_CGATAC192_HCC44_Rep1_S15_R2_001.fastq.gz",
    "20230817_CGATAC193_HCC44_Rep3_S16_R1_001.fastq.gz",
    "20230817_CGATAC193_HCC44_Rep3_S16_R2_001.fastq.gz",
    "20230822_CGATAC194_H358_Rep1_S1_R1_001.fastq.gz",
    "20230822_CGATAC194_H358_Rep1_S1_R2_001.fastq.gz",
    "20230822_CGATAC195_H358_Rep3_S2_R1_001.fastq.gz",
    "20230822_CGATAC195_H358_Rep3_S2_R2_001.fastq.gz",
    "20230822_CGATAC196_H1373_Rep1_S3_R1_001.fastq.gz",
    "20230822_CGATAC196_H1373_Rep1_S3_R2_001.fastq.gz",
    "20230822_CGATAC197_H1373_Rep3_S4_R1_001.fastq.gz",
    "20230822_CGATAC197_H1373_Rep3_S4_R2_001.fastq.gz",
    "20230822_CGATAC198_HOP62_Rep2_S5_R1_001.fastq.gz",
    "20230822_CGATAC198_HOP62_Rep2_S5_R2_001.fastq.gz",
    "20230822_CGATAC199_HOP62_Rep3_S6_R1_001.fastq.gz",
    "20230822_CGATAC199_HOP62_Rep3_S6_R2_001.fastq.gz",
    "20230822_CGATAC200_H23_Rep1_S7_R1_001.fastq.gz",
    "20230822_CGATAC200_H23_Rep1_S7_R2_001.fastq.gz",
    "20230822_CGATAC201_H23_Rep3_S8_R1_001.fastq.gz",
    "20230822_CGATAC201_H23_Rep3_S8_R2_001.fastq.gz",
    "20230822_CGATAC202_Lu65_Rep1_S9_R1_001.fastq.gz",
    "20230822_CGATAC202_Lu65_Rep1_S9_R2_001.fastq.gz",
    "20230822_CGATAC203_Lu65_Rep3_S10_R1_001.fastq.gz",
    "20230822_CGATAC203_Lu65_Rep3_S10_R2_001.fastq.gz",
    "20230822_CGATAC204_H2030_Rep1_S11_R1_001.fastq.gz",
    "20230822_CGATAC204_H2030_Rep1_S11_R2_001.fastq.gz",
    "20230822_CGATAC205_H2030_Rep3_S12_R1_001.fastq.gz",
    "20230822_CGATAC205_H2030_Rep3_S12_R2_001.fastq.gz",
    "20230822_CGATAC206_H2122_Rep1_S13_R1_001.fastq.gz",
    "20230822_CGATAC206_H2122_Rep1_S13_R2_001.fastq.gz",
    "20230822_CGATAC207_H2122_Rep3_S14_R1_001.fastq.gz",
    "20230822_CGATAC207_H2122_Rep3_S14_R2_001.fastq.gz",
    "20230822_CGATAC208_HCC44_Rep1_S15_R1_001.fastq.gz",
    "20230822_CGATAC208_HCC44_Rep1_S15_R2_001.fastq.gz",
    "20230822_CGATAC209_HCC44_Rep3_S16_R1_001.fastq.gz",
    "20230822_CGATAC209_HCC44_Rep3_S16_R2_001.fastq.gz"
]

fastqs = [Path(f) for f in fastqs]

ss1 = pd / 'claudia_sample_sheets' / 'ss1.csv'
ss2 = pd / 'claudia_sample_sheets' / 'ss2.csv'

d1 = get_model_from_sample_sheet(ss1, verbose=False, condition_slice=slice(2, 3))
d2 = get_model_from_sample_sheet(ss2, verbose=False, condition_slice=slice(2, 3))

print(d1)
print(d2)

combined = combine_runs(d1, d2)

print(combined)

combined_merged_reps = combined.merge_reps()

print(combined_merged_reps)

a = combined_merged_reps.get_fastq_groupings(fastqs)
b = combined_merged_reps.get_ordered_conditions() * 2
c = combined_merged_reps.get_ordered_reps() * 2

for fastq, cond, rep in zip(a, b, c):
    samples = [str(f).split('_')[1] for f in fastq]
    read_dirs = [str(f).split('.')[0].split('_')[-2] for f in fastq]
    print(samples, read_dirs, str(cond), str(rep))

# a = combined.get_fastq_groupings(fastqs)
# b = combined.get_ordered_conditions() * 2
# c = combined.get_ordered_reps() * 2
#
# for fastq, cond, rep in zip(a, b, c):
#     print(str(fastq).split('_')[1], str(cond), str(rep))

merge_fastqs_grouped(a, b, c, Path('.'), verbose=True)






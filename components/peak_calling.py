from pathlib import Path
from typing import List

from utils.job_formatter import ExecParams
from utils.job_manager import JobManager
from utils.path_manager import PathManager, cmdify
from wrappers.macs2 import MACSPeakCaller


class PeakCaller:
    """
    Call peaks on a number of input bams. Creates directories to store peaks
    and intermediary bam files.
    """

    jobs: JobManager
    light_job: ExecParams
    heavy_job: ExecParams
    pm: PathManager

    def __init__(self, job_manager: JobManager,
                 path_manager: PathManager,
                 light_job_params: ExecParams,
                 heavy_job_params: ExecParams,
                 lazy: bool = True) -> None:
        self.pm = path_manager
        self.jobs = job_manager
        self.light_job = light_job_params
        self.heavy_job = heavy_job_params
        self.macs2 = MACSPeakCaller()

        if lazy:
            self.run = job_manager.execute_lazy
        else:
            self.run = job_manager.execute

    def run_peak_calling(self, bams: List[Path],
                         peaks_dir_name: str,
                         cutoff_q_val: float,
                         control_bam: Path = None,
                         *args, **kwargs) -> List[Path]:
        """
        Runs peak-calling, creating
        :param bams: The bams on which to run peak calling.
        :param peaks_dir_name: The name of this batch of peak calling.
        :param cutoff_q_val: Q-value threshold for peakcalling.
        :param control_bam: Path to the bam to use as a control.
        :param args: Arguments to be passed to the macs2.call_peaks() function.
        :param kwargs: Arguments to be passed to the macs2.call_peaks()
        function.
        :return: Paths to the created peak files.
        """
        out_dir = self.pm.project_dir / peaks_dir_name

        self.run(cmdify('mkdir', out_dir),
                 self.light_job)

        peakfiles = []
        for bam in bams:
            filename = str(bam).split('/')[-1]
            peaksfile_name = filename[:-4]
            params = self.heavy_job

            self.run(self.macs2.call_peaks(bam,
                                           cutoff_q_val,
                                           out_dir,
                                           params,
                                           control_bam,
                                           outfile_name=peaksfile_name,
                                           *args, **kwargs), params)
            peakfiles.append(out_dir / (peaksfile_name +
                                        '.bed_peaks.narrowPeak'))

        return peakfiles












import os
from pathlib import Path
from typing import Dict

from constants.program_paths import progs
from utils.exceptions import UpstreamPipelineError
from utils.job_formatter import ExecParams, Job
from utils.path_manager import cmdify, PathManager
from wrappers.wrapper import ProgramWrapper


class Demultiplexer(ProgramWrapper):
    """
    === Description ===
    Manages the demultiplexing of Illumina bcl files.
    https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf
    """

    def _get_dependencies(self) -> Dict[str, str]:
        return {}

    def __init__(self, path_manager: PathManager):
        self._pm = path_manager

    def get_default_demultiplex_job(self, exec_params: ExecParams) -> Job:
        """
        Demultiplex using standard file locations.

        Input:
        project_dir/
            sequencing/
                samplesheet.txt
                <YOUR COMPLICATED SEQUENCING RUN NAME>/

        Output:
        project_dir/
            demult/
                Stats/
                Reports/
                <Any number of fastq files>


        :return:
        """

        # Check that sequencing directory exists and contains the required files.
        seq_dir = self._pm.project_dir / 'sequencing'
        if not seq_dir.exists():
            raise UpstreamPipelineError("Sequencing directory not present.")

        # Check that samplesheet exists. need only contain the word "sample".
        sample_inds = ['sample' in s.lower() for s in os.listdir(seq_dir)]
        samplesheet_exists = sum(sample_inds) == 1
        if not samplesheet_exists:
            raise UpstreamPipelineError("Sample sheet not present in directory.")

        # Check that there are only two file in the seq_results directory
        if not len(sample_inds) == 2:
            raise UpstreamPipelineError("Too Many files in sequencing results directory.")

        seq_res_ind = sample_inds.index(False)
        seq_res_path = seq_dir / os.listdir(seq_dir)[seq_res_ind]

        if not seq_res_path.is_dir():
            raise UpstreamPipelineError('Sequencing Results file is not a directory.')

        out_path = self._pm.project_dir / 'demultiplex'
        if not out_path.exists():
            os.mkdir(out_path)

        sample_ind = sample_inds.index(True)
        samplesheet_path = seq_dir / os.listdir(seq_dir)[sample_ind]

        return self.get_demultiplex_cmd(seq_dir, samplesheet_path, out_path, exec_params)

    def get_demultiplex_cmd(self, sequencing_dir: Path,
                            sample_sheet_path: Path, output_dir: Path,
                            num_cores, reports_dir: Path = None,
                            stats_dir: Path = None, io_cores: int = 4
                            ) -> str:
        """
        Demultiplexes the sequencing data at <sequencing_dir>, outputting fastqs, stats, and reports to
        <output_dir>.

        :param io_cores: The number of cores to use for reading and writing.
            cannot (and most certaintly should not) be more than # of
            processing cores. Adding too many of these can lead to segmentation
            faults.
        :param stats_dir: Directory into which the stats should be placed
        :param reports_dir: Directory into which the reports should be placed.
        :param sequencing_dir: Directory containing results of Illumina
            sequencing run.
        :param sample_sheet_path: The path to the sample sheet describing the
            indexing regions used in the given run.
        :param output_dir: The directory into which raw fastqs will be placed,
            alongside run reports and stats.
        :param num_cores: The number of cores to use for processing bcl data.
        :return: A job that will execute the demultiplexing procedure as described.
        """

        io_cores = min([io_cores, num_cores])

        optionals = []
        if reports_dir:
            optionals.extend(['--reports-dir', reports_dir])
        if stats_dir:
            optionals.extend(['--stats-dir', stats_dir])

        cmd = cmdify(
            progs.bcl2fastq,
            '-r', io_cores,
            '-p', num_cores,
            '-w', io_cores,
            '--runfolder-dir', sequencing_dir,
            '--output-dir', output_dir,
            '--sample-sheet', sample_sheet_path,
            "--ignore-missing-bcls",
            "--no-lane-splitting",
            "--barcode-mismatches \"1,1\"",
            *optionals
        )
        return cmd





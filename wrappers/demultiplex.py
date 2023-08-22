import os
from pathlib import Path
from constants.program_paths import progs
from utils.exceptions import UpstreamPipelineError
from utils.job_formatter import ExecParams, Job
from utils.path_manager import cmdify, PathManager


class Demultiplexer:
    """
    === Description ===
    Manages the demultiplexing of Illumina bcl files.
    https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf

    === Private Attributes ===
    """

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

        return self.get_demultiplex_job(seq_dir, samplesheet_path, out_path, exec_params)

    def get_demultiplex_job(self, sequencing_dir: Path, sample_sheet_path: Path, output_dir: Path,
                            exec_params: ExecParams) -> Job:
        """
        Demultiplexes the sequencing data at <sequencing_dir>, outputting fastqs, stats, and reports to
        <output_dir>.

        :param sequencing_dir: Directory containing results of Illumina sequencing run.
        :param sample_sheet_path: The path to the sample sheet describing the indexing regions used in the given run.
        :param output_dir: The directory into which raw fastqs will be placed, alongside run reports and stats.
        :param exec_params: The parameters used to execute the given command.
        :return: A job that will execute the demultiplexing procedure as described.
        """
        num_cores = exec_params.num_cores
        cmd = cmdify(
            progs.bcl2fastq,
            '-r', num_cores,
            '-p', num_cores,
            '-w', num_cores,
            '--runfolder-dir', sequencing_dir,
            '--output_dir', output_dir,
            '--sample_sheet', sample_sheet_path
        )
        return exec_params.builder.prepare_job(cmd, exec_params)





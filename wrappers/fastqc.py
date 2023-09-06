from pathlib import Path
from typing import List, Dict

from utils.job_formatter import ExecParams, Job
from utils.path_manager import cmdify
from wrappers.wrapper import ProgramWrapper


class FastQC(ProgramWrapper):
    """
    === Description ===
    Run quality control on fastq files.
    """

    def _get_dependencies(self) -> Dict[str, str]:
        return {'fastqc': '0.11.9'}

    def fast_qc_arrayed(self, fastq_files: List[Path], out_dir: Path,
                            exec_params: ExecParams) -> List[Job]:
        """
        Runs fastqc, starting enough slurm jobs to assign each fastq file a
        core.
        :param fastq_files: A list of fastq files.
        :param out_dir: The output directory.
        :param exec_params:
        :return:
        """
        n_cores = exec_params.num_cores
        fqc = []
        jobs = []

        for fastq in fastq_files:
            fqc.append(fastq)
            if len(fqc) == n_cores:
                cmd = self.run_fast_qc(fqc, out_dir, exec_params)
                jobs.append(exec_params.builder.prepare_job(cmd, exec_params))
                fqc = []

        if fqc:
            cmd = self.run_fast_qc(fqc, out_dir, exec_params)
            jobs.append(exec_params.builder.prepare_job(cmd, exec_params))

        return jobs


    def run_fast_qc(self, fastq_files: List[Path], out_dir: Path,
                    exec_params: ExecParams) -> str:
        """
        Creates a command that would perform fastqc on the given fasta files.
        :param fastq_files: A list of fastq files.
        :param out_dir: The directory in which to place the resultant fastqc
            results.
        :param threads: The number of threads to use.
        :param exec_params: The parameters required to execute
        :return: A command that, when executed will perform the requested task.
        """

        cmd = cmdify(
            'fastqc',
            *fastq_files,
            '-o', out_dir,
            '-t', exec_params.num_cores
        )

        exec_params.requires.update(self._get_dependencies())

        return cmd






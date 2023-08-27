from abc import ABC
from pathlib import Path
from typing import Dict

from utils.job_formatter import Job, ExecParams
from utils.path_manager import cmdify
from wrappers.wrapper import ProgramWrapper


class MACSPeakCaller(ProgramWrapper):
    """
    === Description ===
    Calls peaks using the MACS2 algorithm.
    """

    requires = {'gcc': '6.2.0',
                'python': '2.7.12',
                'macs2': '2.1.1.20160309'}

    def _get_dependencies(self) -> Dict[str, str]:
        return self.requires

    def call_peaks(self, bam_file: Path,
                   cutoff_q_val: float,
                   output_directory: Path,
                   params: ExecParams,
                   control_bam: Path = None,
                   paired_end: bool = False,
                   genome: str = 'hs',
                   peak_type: str = 'broad',
                   outfile_name: str = None,
                   use_shift_model: bool = True,
                   extension_size: int = 150,
                   broad_join_cutoff: float = 0.05) -> Job:
        """
        Calls peaks on the given bam file, producing a bed file of peaks.
        :param bam_file: The bam file containing the aligned reads of interest.
        :param control_bam: The bam generated from the "input" control, without pulling down for the POI.
        :param paired_end: Whether the given reads are paired-end.
        :param cutoff_q_val: The q value used to decide which peaks are significant.
        :param output_directory: The into which the resulting bed file should be placed.
        :param params: The ExecParams used to configure the Job.
        :param genome: The genome to be used. One of {hs,mm,vero}.
        :param peak_type: Whether to merge nearby peaks into a single broader peak.
        :param outfile_name: The name of the file to write out. Name of the bam file by default.
        :param use_shift_model: Whether to use the shift model, as described in the MACS2 documentation.
        :param extension_size: How far to extend reads when modelling fragments - used in lieu of shift model.
        :param broad_join_cutoff: q-value for how close together peaks have to be to be joined into a single peak.
        :return: Job that will execute the MACS2 script as described above.
        """

        if paired_end:
            pe_str = "BAMPE"
        else:
            pe_str = "BAM"

        optionals = []

        if control_bam:
            optionals.extend(["-c", control_bam,])

        if outfile_name is not None:
            optionals.extend(["-n", outfile_name])

        if not use_shift_model:
            optionals.extend(['--nomodel', "--extsize", extension_size])

        if peak_type == 'broad':
            optionals.append("--broad --broad-cutoff " + str(broad_join_cutoff))

        optionals = tuple(optionals)

        cmd = cmdify("macs2", "callpeak",
                     "-t", bam_file,
                     "-f", pe_str,
                     "-q", cutoff_q_val,
                     "-g", genome,
                     "--outdir", output_directory,
                     *optionals)

        params.requires.update(self.requires)
        return params.builder.prepare_job(cmd, exec_params=params)
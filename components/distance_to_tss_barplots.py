import os
from copy import copy
from pathlib import Path
from typing import Callable, List, Dict

from utils.fetch_files import get_matching_files, get_unique_filename
from utils.job_formatter import ExecParams
from utils.job_manager import JobManager
from utils.path_manager import PathManager, cmdify


class DistanceToTSS:
    """

    === Description ===
    Responsible for augmenting and extracting the distance to TSS stats from
    bed files.

    === Private Attributes ===
    jobs_manager: Handle how compute-intensive jobs are handled.
    path_manager: Manges where this class places files.
    start_array: Command to execute in order to initialise a job array.
    wait_array: Command to execute in order to wait for all jobs in the array
        to complete.
    light_job: Parameters to execute in order to execute a near-instantaneous
        job.
    med_job: Parameters to execute in order to execute a medium length job.


    """
    _jobs: JobManager
    _paths: PathManager
    _start_array: Callable
    _wait_array: Callable
    _light_job: ExecParams
    _med_job: ExecParams

    def __init__(self, jobs_manager: JobManager, path_manager: PathManager,
                 start_array: Callable, wait_array: Callable,
                 light_job: ExecParams, med_job: ExecParams) -> None:
        """
        Initialises this class with the given attributes.
        :param jobs_manager: Handle how compute-intensive jobs are handelled.
        :param path_manager: Manges where this class places files.
        :param start_array: Command to execute in order to initialise a job
            array.
        :param wait_array: Command to execute in order to wait for all jobs in
            the array to complete.
        :param light_job: Parameters to execute in order to execute a
            near-instantaneous job.
        :param med_job: Parameters to execute in order to execute a medium
            length job.
        """
        self._jobs = jobs_manager
        self._paths = path_manager
        self._start_array = start_array
        self._wait_array = wait_array
        self._light_job = copy(light_job)
        self._med_job = copy(med_job)
        self._med_job.wait = True
        self._light_job.wait = True

    def add_positional_info(self, beds_dir: Path, out_dir: Path) -> List[Path]:
        """
        For each bed in beds_dir, creates a txt file that adds the position of
        the nearest gene and the distance to the nearest transcriptional start
        site in out_dir.
        :param beds_dir: The directory in which the beds to be augmented are
            contained.
        :param out_dir: The directory into which the augmented beds should be
            placed.
        :return: Paths to each of the nearest gene beds created.
        """

        beds = get_matching_files(beds_dir, [".*\.bed$"], paths=True)
        out_files = [beds_dir / (str(bed)[:-4] + ".nearestGene.txt")
                     for bed in beds]
        os.chdir(beds_dir.parent)
        cmd = cmdify("perl $soft/addNearestGeneToBED.pl",
                     "$soft/hg19.ensembl.genebody.protein_coding.txt",
                     str(beds_dir) + '/')
        self._jobs.execute_lazy(cmd, self._med_job)
        os.chdir(self._paths.project_dir)

        self._jobs.execute_lazy(cmdify('mv', *out_files, out_dir))

        return [out_dir / out_file.name for out_file in out_files]

    def make_distance_to_tss_tsv(self, nearest_gene_beds: List[Path],
                                 distance_to_tss_tsv: Path,
                                 tmpfile_name: str = None) -> None:
        """
        Creates a distance to TSS tsv containing the distance to TSS stats
        for each of the beds in nearest_gene_beds_dir.
        :param nearest_gene_beds_dir: Directory containing each of the
            augmented nearest gene beds to be included in the tss matrix.
        :param distance_to_tss_tsv: TSV file containing information about
            the distance to the nearest TSS for each bed in
            <nearest_gene_beds_dir>.
        """
        temp_dir = self._paths.purgeable_files_dir
        self._jobs.execute_lazy(cmdify('mkdir', temp_dir))
        self._jobs.execute_lazy(cmdify('cp', *nearest_gene_beds, temp_dir))

        if not tmpfile_name:
            tmpfile_name = distance_to_tss_tsv.name + '.tmp'

        os.chdir(temp_dir.parent)
        cmd = cmdify("perl $soft/getDistanceToTSSmatrix.pl",
                     str(temp_dir) + '/',
                     tmpfile_name,
                     "\n",
                     'mv', tmpfile_name, distance_to_tss_tsv)
        self._jobs.execute_lazy(cmd)
        self._jobs.execute_lazy(cmdify('rm -r', temp_dir))

        os.chdir(self._paths.project_dir)

    def generate_tss_barplot(self, bed_to_bar_name: Dict[Path, str],
                             storage_dir: Path,
                             figure_out_path: Path = None) -> Path:
        """
        Generates stacked TSS bargraphs for the given bedfiles.
        :param beds: Bed files to include in the barplot.
        :param storage_dir: Directory in which to store intermediary files.
        :param bed_to_bar_name: The name for each bedfile to be displayed in
            the barplot.
        :param figure_out_path: Path to final figure. 'out.svg' by default.
        :return: The path to the final figure.
        """
        beds_dir = storage_dir / 'raw_beds'
        positional_info_dir = storage_dir / 'bed_with_pos_info'

        self._jobs.execute_lazy(
            cmdify('mkdir', beds_dir, positional_info_dir))

        for bed in bed_to_bar_name:
            bar_name = bed_to_bar_name[bed]
            self._jobs.execute_lazy(cmdify('cp', bed,
                                           beds_dir / (bar_name + '.bed')))

        aug_beds = self.add_positional_info(beds_dir, positional_info_dir)

        distance_to_tss_tsv = storage_dir / 'distance_to_tss.tsv'
        self.make_distance_to_tss_tsv(aug_beds, distance_to_tss_tsv)

        if not figure_out_path:
            figure_out_path = storage_dir / 'out.svg'
        print(*bed_to_bar_name.values())
        cmd = cmdify('Rscript $Rcode/atac_seq/make_distance_to_tss_barplot.R',
                     distance_to_tss_tsv, figure_out_path,
                     *bed_to_bar_name.values())
        self._jobs.execute_lazy(cmd)

        return figure_out_path





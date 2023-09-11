from copy import copy
from pathlib import Path
from typing import List, Callable

import numpy as np

from utils.fetch_files import get_unique_filename
from utils.job_formatter import ExecParams
from utils.job_manager import JobManager
from utils.path_manager import PathManager, cmdify
from constants.data_paths import HG19_IDXSTATS

HG19_GENOME = HG19_IDXSTATS

class PeakPCAAnalyser:
    """

    === Description ===
    Perform principal component analysis on the peaks of some number of bam and
    bed files.

    === Private Attributes ===
    jobs: The job manager used to execute commands.
    files: The file manager used to find useful paths.
    light_job: Params used to execute near-instantaneous commands.
    parallel_job: Params used to execute long-running jobs that benefit from
        additional threading.
    heavy_job: Params used to execute long-running jobs.

    """

    _jobs: JobManager
    _files: PathManager

    _light_job: ExecParams
    _parallel_job: ExecParams
    _heavy_job: ExecParams

    _start_job_array: Callable
    _wait_job_array: Callable

    def __init__(self, job_manager: JobManager, file_manager: PathManager,
                 light_job_params: ExecParams, heavy_job_params: ExecParams,
                 parallel_job_params: ExecParams,
                 start_job_array_fxn: Callable, wait_job_array_fxn: Callable) \
            -> None:
        """
        Initialises a PCA analyser that executes its jobs using the given
        attributes.

        :param job_manager:
        :param file_manager:
        :param light_job_params: Parameters for a very light job that is
        practically instantaneous.
        :param heavy_job_params: Parameters for a heavy job that does not
            benefit from parallelization. This should have at least 2GB of ram
            and 1 core.
        :param parallel_job_params: Parameters for a heavy job that will
            benefit from parallelization.
        :param start_job_array_fxn: Function to call when starting an array of
            jobs.
        :param wait_job_array_fxn: Function to call in order to wait for an
            array of jobs.
        """

        self._jobs = job_manager
        self._files = file_manager
        self._light_job = copy(light_job_params)
        self._parallel_job = copy(heavy_job_params)
        self._heavy_job = copy(parallel_job_params)
        self._heavy_job.wait = True
        self._heavy_job.add_requirements({"bedtools": "2.30.0"})

        self._start_job_array = start_job_array_fxn
        self._wait_job_array = wait_job_array_fxn

    def find_all_common_peaks(self, beds: list[Path], common_peak_out: Path,
                              genome_index: Path = HG19_GENOME) -> None:
        """
        Merges all the given beds into a single sorted bedfile.

        === Input Data Assumptions ===

        :param beds: A list of bedfiles.
        :param genome_index: A genome index file.
        :param common_peak_out: Path in which to create the merged bed file.
        :return:
        """
        cat_unmerged = self._files.puregable_files_dir / (get_unique_filename() + '.bed')
        self._jobs.execute_lazy(cmdify('cat', *beds, '>', cat_unmerged),
                                self._heavy_job)

        merged = self._files.puregable_files_dir / (
                    get_unique_filename() + '.bed')

        self._jobs.execute_lazy(
            cmdify('bedtools merge -i', cat_unmerged, '>', merged),
            self._heavy_job)

        self._jobs.execute_lazy(
            cmdify('bedtools sort',
                   '-i', merged,
                   '-faidx', genome_index,
                   '>', common_peak_out),
            self._heavy_job)

    def generate_counts(self, bedfile: Path, bamfiles: List[Path],
                        out_dir: Path, out_count_names: List[str] = None,
                        genome_index: Path = HG19_GENOME) -> List[Path]:
        """
        Gets the number of reads in each bamfile in <bamfiles> that lie at each
        peak in bedfile in bedfiles.

        === Input Data Assumptions ===
        The input bedfile and bamfiles must be sorted.
        All input bamfiles have unique characters preceding the first '.'.

        :param genome_index: The genome index for use in when running the
            intersection.
        :param bedfile: A sorted bedfile.
        :param bamfiles: Some number of bam files aligned to the same genome
            as the given bedfile.
        :param out_dir: The directory into which the counts files will be
        placed.
        :return: The paths to the output counts files.
        """

        self._start_job_array()
        new_count_names =[]
        for i, bamfile in enumerate(bamfiles):

            if out_count_names:
                count_file_name = out_count_names[i]
            else:
                count_file_name = bamfile.name.split('.')[0] + '.cnt'

            if count_file_name in new_count_names:
                raise ValueError('Count file collision detected, '
                                 'please check the format of you input bams '
                                 'or explicitly state count file names.')
            self._jobs.execute_lazy(
                cmdify(
                    'bedtools intersect',
                    '-a', bedfile,
                    '-b', bamfile,
                    '-g', genome_index,
                    '-c',
                    '-sorted',
                    '>', out_dir / count_file_name
                ), self._heavy_job
            )
            new_count_names.append(count_file_name)

        self._wait_job_array()

        return [out_dir / count_file_name
                for count_file_name in new_count_names]

    def make_counts_matrix(self, counts_files: List[Path],
                           matrix_out_path: Path,
                           counts_names: List[str] = None) -> np.array:
        """
        Aggreagates all the counts in the given counts files into a single
        matrix

        === Input Data Assumptions ===
        if counts_names is not provided, counts_file.name.split('_')[1] will
        be used as the column name for each counts file.

        counts files have four columns, containing the:
        contig, feature_start, feature_stop, and counts, e.g.

        File 1:
            chr1	954905	955391	12
            chr1	968272	968833	900
            chr1	1002335	1002687	55
            chr1	1051482	1052010	12
            chr1	1167309	1167619	26
            chr1	1209066	1209536	1
            chr1	1243413	1243987	0
            ...

        File 2
            chr1	954905	955391	45
            chr1	968272	968833	429
            chr1	1002335	1002687	44
            chr1	1051482	1052010	58
            chr1	1167309	1167619	26
            chr1	1209066	1209536	51
            chr1	1243413	1243987	53
            ...

        === Output ===

        SAMPLE1, SAMPLE2, SAMPLE 3
        1   4    0
        0   1    1
        1   1    1
        900 50  50
        ...

        Creates a csv file for


        :param counts_files: A list of count bed files, with the counts in the
            fourth column. All counts must have been generated from a single
            bed file.
        :param matrix_out_path: The path to the matrix tsv file to create.
        :param counts_names: The name of each column in the counts_files array.
        :return:
        """

        with open(counts_files[0], 'r') as peaksfile:
            nrow = len(peaksfile.readlines())
        ncol = len(counts_files)

        counts_array = np.zeros((nrow, ncol), dtype=int)

        parsed_col_names = []

        for i, counts_file in enumerate(counts_files):
            if counts_names:
                col_name = counts_names[i]
                parsed_col_names.append(col_name)
            else:
                col_name = counts_file.name.split('_')[1]
                parsed_col_names.append(col_name)

            with open(counts_file, 'r') as counts_file_obj:
                lines = counts_file_obj.readlines()

            for j, vals in enumerate(lines):
                if not vals.strip():
                    continue

                vals = vals.split('\t')

                contig, start, stop, count = vals
                counts_array[j, i] = int(count)

        with open(matrix_out_path, 'w') as outfile:
            outfile.write('\t'.join(parsed_col_names) + '\n')

            for row_num in range(nrow):
                outfile.write('\t'.join(list(counts_array[:, row_num])) + '\n')

        return counts_array

    def generate_pca_plot(self, counts_matrix: Path, *args, **kwargs) -> None:
        """
        Generates a PCA plot displaying a dimensionality reduced -
        representation of the given counts matrix.

        :param counts_matrix: A matrix with the number of counts for each
            peak.
        :param args: Additional argumdents to pass to the graphing function.
        :param kwargs: Additional arguments to pass to the graphing function.
        :return:
        """
        pass


    def do_peaks_pca_analysis(self, beds: List[Path], bams: List[Path],
                              analysis_dir: Path) -> None:
        """
        Performs PCA analysis on the given bedfiles, placing results and
        intermediary files in <analysis_dir>.


        :param beds: A list of bedfiles, to find common peaks on and generate
            PCA plots for.
        :param bams: A list of bam files, to aggregate counts from.
        :param analysis_dir: The directory in which to place results and
            intermediary files.
        """
        if not analysis_dir.exists():
            self._jobs.execute(cmdify('mkdir', analysis_dir))

        common_peaks_path = analysis_dir / 'common_peaks'
        self.find_all_common_peaks(beds, common_peaks_path)

        counts_dir = analysis_dir / 'counts'
        if not counts_dir.exists():
            self._jobs.execute(cmdify('mkdir', counts_dir))

        counts_files = self.generate_counts(common_peaks_path, bams,
                                            counts_dir)

        matrix_path = analysis_dir / 'counts_matrix.tsv'
        self.make_counts_matrix(counts_files, matrix_path)

        self.generate_pca_plot()






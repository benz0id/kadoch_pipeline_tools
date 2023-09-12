import itertools
import logging
import os
import subprocess
from copy import copy
from pathlib import Path
from typing import List, Callable

import numpy as np
import scipy

from utils.cache_manager import CacheManager
from utils.fetch_files import get_unique_filename
from sklearn.preprocessing import StandardScaler
from utils.job_formatter import ExecParams, PythonJob
from utils.job_manager import JobManager
from utils.path_manager import PathManager, cmdify
from constants.data_paths import HG19_IDXSTATS
from utils.utils import ExperimentalDesign

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from utils.utils import ExperimentalDesign

HG19_GENOME = HG19_IDXSTATS

logger = logging.getLogger(__name__)


def generate_pca_plot(counts_matrix_path: Path,
                      design: ExperimentalDesign) -> sns.scatterplot:
    """
    Generates a PCA plot displaying a dimensionality reduced -
    representation of the given counts matrix.

    :param design: The design of the experiment.
    :param counts_matrix_path: A matrix with the number of counts for each
        peak.
    :return:
    """
    # Extract raw data from the counts matrix.
    # Extract raw data from the counts matrix.
    counts_dataframe = pd.read_csv(counts_matrix_path, sep='\t')
    sample_names = counts_dataframe.columns
    counts_matrix = counts_dataframe.loc[:, :].values

    print(counts_dataframe.loc[:40, :].values)


    # Normalise sample.
    # norm_counts = scipy.stats.zscore(counts_matrix.T)
    norm_counts = counts_matrix.T

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(norm_counts)

    pcdf = pd.DataFrame(data=pcs,
                        columns=['principal component 1',
                                 'principal component 2'])

    pcdf['reps'] = [design.get_rep_num(sample)
                    for sample in sample_names]
    pcdf['labels'] = [design.get_condition(sample)
                      for sample in sample_names]

    props = pca.explained_variance_ratio_*100

    ax = sns.scatterplot(data=pcdf,
                         x='principal component 1',
                         y='principal component 2',
                         hue='labels', style='reps',
                         palette=sns.cubehelix_palette(
                             len(design.get_conditions())))

    ax.set(title='PCA',
           xlabel='PC%d:%.2f%%' % (1, props[0]),
           ylabel='PC%d:%.2f%%' % (2, props[1])
           )

    return ax


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

    _verbose: bool

    def __init__(self, job_manager: JobManager, file_manager: PathManager,
                 light_job_params: ExecParams, heavy_job_params: ExecParams,
                 parallel_job_params: ExecParams,
                 start_job_array_fxn: Callable, wait_job_array_fxn: Callable,
                 verbose: bool = False) \
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
        self._parallel_job = copy(parallel_job_params)
        self._heavy_job = copy(heavy_job_params)
        self._heavy_job.wait = True
        self._heavy_job.add_requirements({"bedtools": "2.30.0"})

        self._start_job_array = start_job_array_fxn
        self._wait_job_array = wait_job_array_fxn
        self._verbose = verbose

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
        unique_filename = '.'.join(str(common_peak_out).split('/'))

        cat_unmerged = self._files.puregable_files_dir / (
                    unique_filename + '.cat.bed')
        self._jobs.execute_lazy(cmdify('cat', *beds, '>', cat_unmerged),
                                self._heavy_job)

        merged = self._files.puregable_files_dir / (
                unique_filename + '.merged.bed')

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
        new_count_names = []
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

    def get_number_of_mapped_reads(self, sample_names: List[str],
                                   bamfiles: List[Path]) \
            -> List[int]:
        """
        Gets the reads that mapped in for each sample in <sample names>
        For each sample name, there must exist a bam file containing that name
        in <bamfiles>.

        :param sample_names: The names of a given number of samples.
        :param bamfiles: The bamfiles generated for each of those samples
        :return: The number of reads in each bamfile, ordered by the
            <sample_names>.
        """

        # Construct map from each count file to the bam file with a
        # matching name, ensuring a 1:1 mapping.
        count_to_bam_map = {}
        for col_name in sample_names:
            matching_bam = None

            for bam in bamfiles:
                if col_name in bam.name and matching_bam:
                    raise ValueError(
                        "Detected multiple bams with sample names matching"
                        " a count file name.\n\t" + '\n\t'.join([
                            str(col_name), bam.name, matching_bam.name
                        ]))
                elif col_name in bam.name:
                    matching_bam = bam

            if not matching_bam:
                raise ValueError(
                    "The given bam files do not contain a matching name "
                    "for " + str(col_name))
            else:
                count_to_bam_map[col_name] = matching_bam

        norm_factors = []
        for i, col_name in enumerate(sample_names):
            bamfile = count_to_bam_map[col_name]
            result = subprocess.run(['samtools', 'view', '-c', str(bamfile)],
                                    stdout=subprocess.PIPE)
            n_reads = int((result.stdout).decode('utf-8').split(' ')[0])
            norm_factors.append(n_reads)

        s = 'Normalisation info:\n'
        for counts_file, norm_factor in zip(sample_names, norm_factors):
            bam_file = count_to_bam_map[counts_file]
            s += ('\t' + counts_file + '\t' +
                  bam_file.name + '\t' +
                  str(norm_factor) + '\n')

        logger.info(s)
        if self._verbose:
            print(s)

        return norm_factors

    def make_counts_matrix(self, counts_files: List[Path],
                           matrix_out_path: Path,
                           counts_names: List[str] = None,
                           bams_to_normalise_to: List[
                               Path] = None) -> np.array:
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
        12  45    0
        0   1    1
        1   1    1
        900 50  50
        ...



        :param counts_files: A list of count bed files, with the counts in the
            fourth column. All counts must have been generated from a single
            bed file.
        :param matrix_out_path: The path to the matrix tsv file to create.
        :param counts_names: The name of each column in the counts_files array.
        :param bams_to_normalise_to: A list of bamfiles with names matching
            the names of each counts file or the supplied counts names.
        :return:
        """

        with open(counts_files[0], 'r') as peaksfile:
            nrow = len(peaksfile.readlines())
        ncol = len(counts_files)

        counts_array = np.zeros((nrow, ncol), dtype=float)

        parsed_col_names = []

        # Collect all counts values and store in NP array.#TODO imp with pandas
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

        # Normalise to cpms if bam files are provided.
        if bams_to_normalise_to:
            norm_factors = self.get_number_of_mapped_reads(parsed_col_names,
                                                           bams_to_normalise_to)
            for col in range(len(parsed_col_names)):
                fac = norm_factors[col]
                column = counts_array[:, col]
                counts_array[:, col] = column / fac * 10 ** 6

        with open(matrix_out_path, 'w') as outfile:
            outfile.write('\t'.join(parsed_col_names) + '\n')

            for row_num in range(nrow):
                outfile.write('\t'.join(
                    [str(val) for val in counts_array[row_num, :]]) + '\n')

        return counts_array

    def do_peaks_pca_analysis(self, beds: List[Path], bams: List[Path],
                              analysis_dir: Path,
                              experimental_design: ExperimentalDesign = None,
                              normalise_by_counts: bool = True) -> None:
        """
        Performs PCA analysis on the given bedfiles, placing results and
        intermediary files in <analysis_dir>.


        :param experimental_design: The design of this experiment.
        :param normalise_by_counts:
        :param beds: A list of bedfiles, to find common peaks on and generate
            PCA plots for.
        :param bams: A list of bam files, to aggregate counts from.
        :param analysis_dir: The directory in which to place results and
            intermediary files.
        """
        if not analysis_dir.exists():
            self._jobs.execute(cmdify('mkdir', analysis_dir))

        common_peaks_path = analysis_dir / 'common_peaks.bed'
        self.find_all_common_peaks(beds, common_peaks_path)

        counts_dir = analysis_dir / 'counts'
        if not counts_dir.exists():
            self._jobs.execute(cmdify('mkdir', counts_dir))

        counts_files = self.generate_counts(common_peaks_path, bams,
                                            counts_dir)
        counts_files = sorted(counts_files,
                              key=lambda x: str(x).split('_')[1])

        if not normalise_by_counts:
            matrix_path = analysis_dir / 'merged_counts_matrix.tsv'
            bams_to_normalise_to = []
        else:
            matrix_path = analysis_dir / 'merged_counts_matrix_normalized.tsv'
            bams_to_normalise_to = bams

        pj = PythonJob('make ' + str(matrix_path), [], self.make_counts_matrix,
                       counts_files=counts_files,
                       matrix_out_path=matrix_path,
                       bams_to_normalise_to=bams_to_normalise_to)
        self._jobs.execute_lazy(pj)

        plot = generate_pca_plot(matrix_path, experimental_design)
        plt.savefig(analysis_dir / (get_unique_filename() + 'pca.svg'),
                    format='svg')




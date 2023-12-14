import itertools
import logging
import os
import subprocess
from copy import copy
from pathlib import Path
from typing import List, Callable, Tuple

import numpy as np
import qnorm as qnorm
from pandas import unique
from scipy.stats import stats

from utils.cache_manager import CacheManager
from utils.fetch_files import get_unique_filename, outpath_to_dirname
from sklearn.preprocessing import StandardScaler
from utils.job_formatter import ExecParams, PythonJob
from utils.job_manager import JobManager
from utils.path_manager import PathManager, cmdify
from constants.data_paths import HG19_IDXSTATS
from utils.utils import ExperimentalDesign
from operator import itemgetter

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from utils.utils import ExperimentalDesign

HG19_GENOME = HG19_IDXSTATS

logger = logging.getLogger(__name__)

MAX_PCA_DIMS = 5


def generate_pca_plot(counts_matrix_path: Path,
                      design: ExperimentalDesign,
                      out_filepath: Path, dims: int = 3,
                      n_info_cols: int = 0, sample_ids: bool = False,
                      samples: List[str] = None,
                      colour_groups: List[str] = 'by_condition',
                      shape_groups: List[str] = 'by_rep',
                      title: str = 'PCA',
                      verbose: bool = False) -> sns.scatterplot:
    """
    Generates a PCA plot displaying a dimensionality reduced -
    representation of the given counts matrix.

    :param samples:
    :param shape_groups:
    :param colour_groups:
    :param dims:
    :param out_filepath:
    :param design: The design of the experiment.
    :param counts_matrix_path: A matrix with the number of counts for each
        peak.
    :param n_info_cols: The number of leading info columns to ignore.
    :param sample_ids: Whether to expect the column headers to be sample ids
        rather than sample names.
    :return:
    """
    # Extract raw data from the counts matrix.
    counts_dataframe = pd.read_csv(counts_matrix_path, sep='\t')

    matrix_samples = counts_dataframe.columns
    if n_info_cols > 0:
        to_rem = [matrix_samples[i] for i in range(n_info_cols)]
        counts_dataframe = counts_dataframe.drop(columns=to_rem)
        matrix_samples = counts_dataframe.columns
        counts_dataframe = counts_dataframe[matrix_samples].astype(float)

    if sample_ids:
        matrix_samples = [sample.strip().split('_')[1] for sample in matrix_samples]
        counts_dataframe.columns = matrix_samples

    if not samples:
        reps = [design.get_rep_num(label) for label in matrix_samples]
        conds = [design.get_condition(label) for label in matrix_samples]
        samples = matrix_samples
    else:
        reps = [design.get_rep_num(label) for label in samples]
        conds = [design.get_condition(label) for label in samples]
        counts_dataframe = counts_dataframe[samples]

    if isinstance(colour_groups, list):
        assert len(colour_groups) == len(samples)
    if isinstance(shape_groups, list):
        assert len(shape_groups) == len(samples)

    counts_matrix = np.log2(counts_dataframe + 1)

    # norm_counts = scaler.fit_transform(counts_matrix)
    norm_counts = stats.zscore(counts_matrix, axis=0)

    pcs = qnorm.quantile_normalize(norm_counts.T)

    pca = PCA(min(MAX_PCA_DIMS, len(samples)))
    pcs = pca.fit_transform(pcs)

    pcdf = pd.DataFrame(data=pcs,
                        columns=['principal component ' + str(i)
                                 for i in range(1, pcs.shape[1] + 1)])
    pcdf['samples'] = samples

    if colour_groups == 'by_condition':
        pcdf['labels'] = conds
    if shape_groups == 'by_rep':
        pcdf['reps'] = reps

    if isinstance(colour_groups, list):
        pcdf['labels'] = colour_groups
    if isinstance(shape_groups, list):
        pcdf['reps'] = shape_groups

    if verbose:
        print(pcdf.iloc[:, -3:])

    logger.info((f'PCA Dataframe: for {counts_matrix_path}:' + '\n' +
                 str(pcdf)))

    props = pca.explained_variance_ratio_ * 100

    # Compare multiple principle components.
    for i, j in itertools.combinations(range(dims), 2):

        kwargs = {
            "data": pcdf,
            "x": 'principal component ' + str(i + 1),
            "y": 'principal component ' + str(j + 1)
        }

        if colour_groups is not None:
            kwargs['hue'] = 'labels'
            kwargs['palette'] = \
                sns.color_palette('bright', len(unique(pcdf['labels'])))
        if shape_groups is not None:
            kwargs["style"] = 'reps'

        ax = sns.scatterplot(**kwargs)

        ax.set(title=title,
               xlabel=f'PC{i + 1}:{props[i]: .2f}',
               ylabel=f'PC{j + 1}:{props[j]: .2f}'
               )
        filename = ''.join(['pc', str(i + 1), '_vs_', 'pc',
                            str(j + 1) + '.svg'])
        plt.show()
        plt.savefig(out_filepath / filename)
        plt.close()


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

        self._idx_stats = HG19_GENOME
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

    def update_genome_index(self, bam: Path) -> None:
        """
        Creates a new genome from the given bam which will be used for all
        subsequent sorting operations.
        :return:
        """
        idx_stats = self._files.project_dir / (bam.name[:-4] + '.genome')
        self._jobs.execute_lazy(cmdify(
            'samtools idxstats', bam,
            "| awk -v OFS='\t' {'print $1,$2'}",
            "| sed '$ d'"
            '>', idx_stats
        ), self._light_job)

        self._idx_stats = idx_stats

    def find_all_common_peaks(self, beds: List[Path], common_peak_out: Path) \
            -> None:
        """
        Merges all the given beds into a single sorted bedfile.

        === Input Data Assumptions ===

        :param beds: A list of bedfiles.
        :param genome_index: A genome index file.
        :param common_peak_out: Path in which to create the merged bed file.
        :return:
        """

        self._jobs.execute_lazy(
            cmdify('cat', *beds, '|',
                   'bedtools sort',
                   '-faidx', self._idx_stats, '|',
                   'bedtools merge',
                   '>', common_peak_out
                   ))

    def generate_counts(self, bedfile: Path, read_files: List[Path],
                        out_dir: Path, out_count_names: List[str] = None) \
            -> List[Path]:
        """
        Gets the number of reads in each read_files in <read_files> that lie at each
        peak in bedfile in <bedfiles>.

        === Input Data Assumptions ===
        The input bedfile and read files must be sorted.
        All input read files have unique characters preceding the first '.'.

        :param bedfile: A sorted bedfile.
        :param read_files: Some number of bam or bed files aligned to the same
            genome as the given bedfile.
        :param out_dir: The directory into which the counts files will be
        placed.
        :return: The paths to the output counts files.
        """

        self._start_job_array()
        new_count_names = []
        for i, read_file in enumerate(read_files):

            if out_count_names:
                count_file_name = out_count_names[i]
            else:
                count_file_name = read_file.name.split('.')[0] + '.cnt'

            if count_file_name in new_count_names:
                raise ValueError('Count file collision detected, '
                                 'please check the format of you input read files '
                                 'or explicitly state count file names.')

            self._jobs.execute_lazy(
                cmdify(
                    'bedtools intersect',
                    '-a', bedfile,
                    '-b', read_file,
                    '-g', self._idx_stats,
                    '-c',
                    '-sorted',
                    '>', out_dir / count_file_name
                ), self._heavy_job
            )
            new_count_names.append(count_file_name)

        self._wait_job_array()

        return [out_dir / count_file_name
                for count_file_name in new_count_names]

    def bam_to_bed(self, bams: List[Path], out_dir: Path,
                   filter_pe: bool = True, sort: bool = True,
                   threads: int = 1, mem: int = 1048,
                   paired_end: bool = True) -> \
            Tuple[List[str], List[Path]]:
        """
        Convert the given bam files into bed files. Place these files into
        <out dir>.


        Sorting will require significant ram and cpu usage.

        :param bams: A list of valid bam files.
        :param out_dir:
        :param filter_pe: Remove improperly paired reads.
        :param sort: Whether to sort bam files before converting them to bed
        files.
        :param threads: The number of threads to spawn when parallelizing.
        :param mem: Max memory usage to use, in MB.
        :param paired_end: Whether to expect paired end input.
        :return: Commands to create beds, paths to created beds.
        """

        old_idx_stats = self._idx_stats
        self.update_genome_index(bams[0])

        out_beds = []
        cmds = []
        sort_str = ''
        filter_arg = ''
        pe_str = ''
        if sort:
            # Half of total memory, divided among threads.
            mem_per = int(mem / 2 / threads)
            sort_str = f'| samtools sort -@ {threads} -m {mem_per:d}M -nu'
        if filter_pe:
            filter_arg = '-f 0x2'
        if paired_end:
            pe_str = '-bedpe'

        if filter_pe and not paired_end:
            raise ValueError('Requested unnecessary paired-end filter.')

        sort_mem = int(mem * 3 / 4)

        for bam in bams:
            out_bed_path = out_dir / (bam.name[:-4] + '.bed')
            tmp_pe_bed = self._files.purgeable_files_dir / \
                         (outpath_to_dirname(out_bed_path)[:-4] + '.bedpe')
            tmp_unsorted_bed = self._files.purgeable_files_dir / \
            (outpath_to_dirname(out_bed_path)[:-4] + '.unsorted.bed')

            cmd = cmdify(
                'samtools view -bu', filter_arg, bam,
                sort_str,
                '| bedtools bamtobed', pe_str, '-i stdin',
                '>', tmp_pe_bed, '\n'
                'cat', tmp_pe_bed,
                "| awk -v OFS='\t' {'print $1,$2,$6,$7,$8,$9'}",
                '>', tmp_unsorted_bed, '\n'
                'bedtools sort -i ', tmp_unsorted_bed, '-faidx', self._idx_stats,
                '>', out_bed_path
            )

            # Old sortbed. Less likely to fail due to memory issues, but does
            # not use faidx.
            # 'sort-bed --max-mem', f'{sort_mem}M', tmp_unsorted_bed,

            out_beds.append(out_bed_path)
            cmds.append(cmd)

        self._idx_stats = old_idx_stats

        return cmds, out_beds

    def get_number_of_mapped_reads(self, sample_names: List[str],
                                   read_files: List[Path]) \
            -> List[int]:
        """
        Gets the reads that mapped in for each sample in <sample names>
        For each sample name, there must exist a bam file containing that name
        in <bamfiles>.

        :param sample_names: The names of a given number of samples.
        :param read_files: The read_files generated for each of those samples
        :return: The number of reads in each read_files, ordered by the
            <sample_names>.
        """

        ft = read_files[0].name.split('.')[-1]
        for file in read_files:
            if not file.name.split('.')[-1] == ft:
                raise ValueError('All read files must be of the same type.')

        # Construct map from each count file to the bam file with a
        # matching name, ensuring a 1:1 mapping.
        count_to_read_map = {}
        for col_name in sample_names:
            matching_read = None

            for read in read_files:
                if col_name in read.name and matching_read:
                    raise ValueError(
                        "Detected multiple reads with sample names matching"
                        " a count file name.\n\t" + '\n\t'.join([
                            str(col_name), read.name, matching_read.name
                        ]))
                elif col_name in read.name:
                    matching_read = read

            if not matching_read:
                raise ValueError(
                    "The given read files do not contain a matching name "
                    "for " + str(col_name))
            else:
                count_to_read_map[col_name] = matching_read

        if ft == 'bam':
            cmd = 'samtools', 'view', '-c'
        elif ft == 'bed':
            cmd = 'wc', '-l'
        else:
            raise ValueError("Unrecognised filetype.")

        norm_factors = []
        for i, col_name in enumerate(sample_names):
            readfile = count_to_read_map[col_name]
            logger.info('Calculating library size for ' + readfile.name)
            result = subprocess.run([*cmd, str(readfile)],
                                    stdout=subprocess.PIPE)
            n_reads = int(result.stdout.decode('utf-8').split(' ')[0])
            norm_factors.append(n_reads)

        s = 'Normalisation info:\n'
        for counts_file, norm_factor in zip(sample_names, norm_factors):
            read_file = count_to_read_map[counts_file]
            s += ('\t' + counts_file + '\t' +
                  read_file.name + '\t' +
                  str(norm_factor) + '\n')

        logger.info(s)
        if self._verbose:
            print(s)

        return norm_factors

    def make_counts_matrix(self, counts_files: List[Path],
                           matrix_out_path: Path,
                           counts_names: List[str] = None,
                           reads_to_normalise_to: List[Path] = None,
                           add_sites_col: bool = False,
                           normalise_by_site_len: bool = False) -> np.array:
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
        12  45   0
        0   1    1
        1   1    1
        900 50  50
        ...

        :param counts_files: A list of count bed files, with the counts in the
            fourth column. All counts must have been generated from a single
            bed file.
        :param matrix_out_path: The path to the matrix tsv file to create.
        :param counts_names: The name of each column in the counts_files array.
        :param reads_to_normalise_to: A list of bamfiles with names matching
            the names of each counts file or the supplied counts names.
        :return:
        """

        with open(counts_files[0], 'r') as peaksfile:
            nrow = len(peaksfile.readlines())
        ncol = len(counts_files)

        counts_array = np.zeros((nrow, ncol), dtype=float)
        sites = [''] * nrow
        site_lens = np.zeros(nrow, dtype=float)

        parsed_col_names = []

        failed = False

        # Collect all counts values and store in NP array.#TODO imp with pandas
        for i, counts_file in enumerate(counts_files):
            logger.info("Starting to Parse Counts from " + str(counts_file))

            if counts_names:
                col_name = counts_names[i]
                parsed_col_names.append(col_name)
            else:
                col_name = counts_file.name.split('_')[1]
                parsed_col_names.append(col_name)

            with open(counts_file, 'r') as counts_file_obj:
                lines = counts_file_obj.readlines()

            for j, vals in enumerate(lines):

                if j % 10000 == 0:
                    logger.debug('\t' + str(j) + ' lines parsed')

                if not vals.strip():
                    continue

                if j >= nrow:
                    failed = True
                    continue

                vals = vals.split('\t')

                contig, start, stop, count = vals
                counts_array[j, i] = int(count)

                site = f'{contig}:{start}-{stop}'

                if not sites[j]:
                    sites[j] = site
                    site_lens[j] = abs(int(stop) - int(start))

        # Normalize to RPMs if bam files are provided.
        if reads_to_normalise_to:
            norm_factors = self.get_number_of_mapped_reads(parsed_col_names,
                                                           reads_to_normalise_to)
            for col in range(len(parsed_col_names)):
                fac = norm_factors[col]
                column = counts_array[:, col]
                counts_array[:, col] = column / fac * 10 ** 6

        # Normalize to RPKMs.
        if normalise_by_site_len:
            site_len_factors = 1 / site_lens * 1000
            counts_array = (counts_array.T / site_len_factors).T

        # Write out matrix
        with open(matrix_out_path, 'w') as outfile:

            if add_sites_col:
                parsed_col_names.insert(0, 'Sites')
            outfile.write('\t'.join(parsed_col_names) + '\n')

            for row_num in range(nrow):
                row = [str(val) for val in counts_array[row_num, :]]
                if add_sites_col:
                    row.insert(0, sites[row_num])

                outfile.write('\t'.join(row) + '\n')
        if failed:
            raise ValueError("Counts files of invalid lengths.")
        return counts_array

    def filter_counts_matrix(self, in_matrix: Path, out_matrix: Path,
                             threshold: int) -> None:
        pass

    def do_peaks_pca_analysis(self, beds: List[Path], reads: List[Path],
                              analysis_dir: Path,
                              experimental_design: ExperimentalDesign = None,
                              normalise_by_counts: bool = True) -> Path:
        """
        Performs PCA analysis on the given bedfiles, placing results and
        intermediary files in <analysis_dir>.


        :param experimental_design: The design of this experiment.
        :param normalise_by_counts: Whether to normalise by library size.
        :param beds: A list of bedfiles, to find common peaks on and generate
            PCA plots for.
        :param reads: A list of read files, either beds or bams, to aggregate
            counts from. Must all be sorted by position using the same genome.
        :param analysis_dir: The directory in which to place results and
            intermediary files.
        """

        ft = reads[0].name.split('.')[-1]
        for file in reads:
            if not file.name.split('.')[-1] == ft:
                raise ValueError('All read files must be of the same type.')

        if ft == 'bam':
            old_idx = self._idx_stats
            self.update_genome_index(reads[0])

        if not analysis_dir.exists():
            self._jobs.execute(cmdify('mkdir', analysis_dir))

        common_peaks_path = analysis_dir / 'common_peaks.bed'
        self.find_all_common_peaks(beds, common_peaks_path)

        counts_dir = analysis_dir / 'counts'
        if not counts_dir.exists():
            self._jobs.execute(cmdify('mkdir', counts_dir))

        counts_files = self.generate_counts(common_peaks_path, reads,
                                            counts_dir)

        def get_pos(f) -> int:
            sample = str(f.name).split('_')[1]
            assert sample in experimental_design.get_samples()
            return experimental_design.get_samples().index(sample)

        counts_files = sorted(counts_files, key=get_pos)

        if not normalise_by_counts:
            matrix_path = analysis_dir / 'merged_counts_matrix.tsv'
            reads_to_normalise_to = []
        else:
            matrix_path = analysis_dir / 'merged_counts_matrix_normalized.tsv'
            reads_to_normalise_to = reads

        counts_files = experimental_design.align_to_samples(counts_files)

        pj = PythonJob('make ' + str(matrix_path), [], self.make_counts_matrix,
                       counts_files=counts_files,
                       matrix_out_path=matrix_path,
                       reads_to_normalise_to=reads_to_normalise_to,
                       add_sites_col=True,
                       normalise_by_site_len=True)
        self._jobs.execute_lazy(pj)

        args = (matrix_path, experimental_design, analysis_dir, 4)
        pj = PythonJob('Generate PCA figures' + str(args[-2:]), [],
                       generate_pca_plot, *args, n_info_cols=1)
        self._jobs.execute_lazy(pj)
        if ft == 'bam':
            self._idx_stats = old_idx

        return matrix_path




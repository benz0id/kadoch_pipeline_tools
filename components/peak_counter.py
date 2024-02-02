import collections
import logging
import subprocess
from copy import copy
from pathlib import Path
from typing import List, Callable, Tuple, Union, Dict

import numpy as np
import threading

from utils.fetch_files import get_unique_filename, outpath_to_dirname, \
    get_matching_files
from utils.job_formatter import ExecParams, PythonJob
from utils.job_manager import JobManager
from utils.path_manager import PathManager, cmdify
from constants.data_paths import HG19_IDXSTATS
from utils.utils import ExperimentalDesign

HG19_GENOME = HG19_IDXSTATS

logger = logging.getLogger(__name__)

file_lock = threading.Lock()

CacheEntry = collections.namedtuple('CacheEntry',
                                    ['bed_file', 'reads_file', 'counts_file'])


class PeakCounter:
    """
    Responsible for generating counts matrices from peak files and reads.
    """

    _jobs: JobManager
    _files: PathManager

    _light_job: ExecParams
    _parallel_job: ExecParams
    _heavy_job: ExecParams

    _start_job_array: Callable
    _wait_job_array: Callable

    _verbose: bool

    _bed_counts_dict: Dict[Path, Path]
    _cache_record_path: Path
    _cache_dir_path: Path

    def __init__(self, job_manager: JobManager, file_manager: PathManager,
                 light_job_params: ExecParams, heavy_job_params: ExecParams,
                 parallel_job_params: ExecParams,
                 start_job_array_fxn: Callable, wait_job_array_fxn: Callable,
                 verbose: bool = False) \
            -> None:
        """
        :param job_manager: Used to execute jobs.
        :param file_manager: Used to navigate the project filesystem.
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
        self._cache_record_path = self._files.cache_dir / 'counts_cache.txt'
        self._cache_dir_path = self._files.make(self._files.cache_dir /
                                           'counts_cache')

    def parse_cache(self, lines: List[str]) -> List[CacheEntry]:
        """
        Parse a cache of reads files, formatted <original_bedfile>\t<counts_bedfile>
        :param lines: A list of lines of the original format.
        :return:
        """
        entries = []
        for line in lines:
            line = line.strip()
            if line == '':
                continue

            components = line.split('\t')
            components = tuple(Path(comp) for comp in components)

            if not len(components) == 3:
                raise RuntimeError(str(self._cache_record_path) +
                                   ' is not properly formatted.')
            entries.append(CacheEntry(*components))
        return entries

    def add_counts_file_to_cache(self, bed: Path, reads_file: Path,
                                 counts: Path) -> None:
        """
        Checks the cache of file for which counts
        :param bed:
        :param reads_file:
        :param counts:
        :return:
        """
        if not bed.exists():
            raise ValueError(str(bed) + ' does not exist.')

        if not reads_file.exists():
            raise ValueError(str(reads_file) + ' does not exist.')

        with file_lock:
            if self._cache_record_path.exists():
                with open(self._cache_record_path, 'r') as cache_file:
                    lines = cache_file.readlines()
                    cache = self.parse_cache(lines)
            else:
                cache = []

            cache.append(CacheEntry(bed, reads_file, counts))

            with open(self._cache_record_path, 'w') as cache_file:
                for entry in cache:
                    cache_file.write(
                        f'{entry.bed_file}\t'
                        f'{entry.reads_file}\t'
                        f'{entry.counts_file}\n')

    def check_cache_for_counts_file(self, bed_file: Path, reads_file: Path) \
            -> Union[None, Path]:
        """
        Checks the cache paths to counts files that have been parsed before.
        :param bed_file: Bed file used to generate counts.
        :param reads_file: Read file used to generate counts.
        :return:
        """
        if not self._cache_record_path.exists():
            return None

        with file_lock, open(self._cache_record_path, 'r') as cache_file:
            lines = cache_file.readlines()
            cache = self.parse_cache(lines)
            for entry in cache:
                bed_match = entry.bed_file == bed_file
                read_match = entry.reads_file == reads_file
                if bed_match and read_match:
                    return entry.counts_file
        return None

    def get_matrix(self, bed: Path, read_file_dir: Path,
                   design: ExperimentalDesign,
                   out_matrix_path: Path, samples: List[str] = None,
                   norm_method: str = 'RPKM') -> None:
        """
        Count the number of reads/fragments at each peak in <bed> for each of
        <read_files>.

        NOTE:
            Columns correspond to ALPHABETICALLY sorted <read_files> names.

        :param bed: A valid bed file. Each annotation in this file will form a
        row in the output matrix.
        :param read_file_dir: A directory containing bed or bam files for
        each of the requested samples.
        :param out_matrix_path: The path into which the parse matrix will be
        placed.
        :param design: Design of the experiment.
        :param samples: The samples to include in the matrix. All samples in
            the provided design by default.
        :param norm_method: How to normalize the reads. Only RPKM is currently
        supported.

        :return:
        """
        if not samples:
            samples = design.get_samples()

        read_files = get_matching_files(
            read_file_dir, samples, filetype='[bam|bed]',
            one_to_one=True, under_delim=True, paths=True)

        # Check that all input files are of the same type.
        ft = read_files[0].name.split('.')[-1]
        for file in read_files:
            if not file.name.split('.')[-1] == ft:
                raise ValueError('All read files must be of the same type.')

        # Update genome index if the file type is a bam.
        if ft == 'bam':
            self.update_genome_index(read_files[0])

        # Sort count files by the experimental design if one was provided.
        column_names = [rf.name.split('_')[1] for rf in read_files]

        # Generate counts files.
        counts_files = self.generate_counts(bed, read_files)

        # Define behaviour based on normalization method.
        if norm_method == 'RPKM':
            add_sites_col = True
            normalise_by_site_len = True
        else:
            raise ValueError('Unrecognised normalization method.')

        # Make the matrix.
        pj = PythonJob('make ' + str(out_matrix_path), [],
                       self.make_counts_matrix,

                       counts_files=counts_files,
                       matrix_out_path=out_matrix_path,
                       reads_to_normalise_to=read_files,
                       column_names=column_names,
                       add_sites_col=add_sites_col,
                       normalise_by_site_len=normalise_by_site_len
                       )
        self._jobs.execute_lazy(pj)

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
                           column_names: List[str] = None,
                           add_sites_col: bool = True,
                           normalise_by_site_len: bool = True,
                           reads_to_normalise_to: List[Path] = None) -> np.array:
        """
        Aggregates all the counts in the given counts files into a single
        matrix. Use of *get_matrix* is generally reccomended for most users.

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

            if column_names:
                col_name = column_names[i]
                parsed_col_names.append(col_name)
            else:
                col_name = counts_file.name
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

    def generate_counts(self, bedfile: Path, read_files: List[Path]) \
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
        :return: The paths to the output counts files.
        """

        self._start_job_array()
        count_files = []

        for i, read_file in enumerate(read_files):

            cached_counts_file = self.check_cache_for_counts_file(bedfile, read_file)
            if cached_counts_file is not None:
                count_files.append(cached_counts_file)
                continue

            out_count_file = self._cache_dir_path / (get_unique_filename() + '.bed')

            self._jobs.execute_lazy(
                cmdify(
                    'bedtools sort',
                    '-i', bedfile,
                    '-g', self._idx_stats,
                    ' | bedtools intersect',
                    '-a stdin',
                    '-b', read_file,
                    '-g', self._idx_stats,
                    '-c',
                    '-sorted',
                    '>', out_count_file
                ), self._heavy_job
            )
            count_files.append(out_count_file)
            self.add_counts_file_to_cache(bedfile, read_file, out_count_file)

        self._wait_job_array()

        return count_files

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
                               (outpath_to_dirname(out_bed_path)[
                                :-4] + '.unsorted.bed')

            cmd = cmdify(
                'samtools view -bu', filter_arg, bam,
                sort_str,
                '| bedtools bamtobed', pe_str, '-i stdin',
                '>', tmp_pe_bed, '\n'
                                 'cat', tmp_pe_bed,
                "| awk -v OFS='\\t' {'print $1,$2,$6,$7,$8,$9'}",
                '>', tmp_unsorted_bed, '\n'
                                       'bedtools sort -i ', tmp_unsorted_bed,
                '-faidx', self._idx_stats,
                '>', out_bed_path
            )

            # Old sortbed. Less likely to fail due to memory issues, but does
            # not use faidx.
            # 'sort-bed --max-mem', f'{sort_mem}M', tmp_unsorted_bed,

            out_beds.append(out_bed_path)
            cmds.append(cmd)

        self._idx_stats = old_idx_stats

        return cmds, out_beds

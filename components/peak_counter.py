
import logging
import subprocess
from copy import copy
from pathlib import Path
from typing import List, Callable, Tuple

import numpy as np

from utils.fetch_files import get_unique_filename, outpath_to_dirname
from utils.job_formatter import ExecParams, PythonJob
from utils.job_manager import JobManager
from utils.path_manager import PathManager, cmdify
from constants.data_paths import HG19_IDXSTATS
from utils.utils import ExperimentalDesign

HG19_GENOME = HG19_IDXSTATS

logger = logging.getLogger(__name__)


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

    def get_matrix(self, bed: Path,  read_files: List[Path],
                   out_matrix_path: Path, norm_method: str = 'RPKM',
                   design: ExperimentalDesign = None) -> None:
        """
        Count the number of reads/fragments at each peak in <bed> for each of
        <read_files>.

        NOTE:
            Columns correspond to ALPHABETICALLY sorted <read_files> names.

        :param bed: A valid bed file. Each annotation in this file will form a
        row in the output matrix.
        :param read_files: A number of BAM/BED files containing read/fragment
        information respectively.
        :param out_matrix_path: The path into which the parse matrix will be
        placed.
        :param norm_method: How to normalize the reads. Only RPKM is currently
        supported.
        :param design: Sort files by their sample name if a design is provided.

        :return:
        """
        # Check that all input files are of the same type.
        ft = read_files[0].name.split('.')[-1]
        for file in read_files:
            if not file.name.split('.')[-1] == ft:
                raise ValueError('All read files must be of the same type.')

        # Update genome index if the filtype is a bam.
        if ft == 'bam':
            self.update_genome_index(read_files[0])

        # Create temporary directory.
        tmp_dir = self._files.make(
            self._files.purgeable_files_dir /
            outpath_to_dirname(out_matrix_path))

        counts_dir = tmp_dir / 'counts'
        self._files.make(counts_dir)

        # Generate counts files.
        counts_files = self.generate_counts(bed, read_files, counts_dir)

        # Define behaviour based on normalization method.
        if norm_method == 'RPKM':
            add_sites_col = True
            normalise_by_site_len = True
        else:
            raise ValueError('Unrecognised normalization method.')

        # Sort count files by the experimental design if one was provided
        if design is not None:
            counts_files = design.align_to_samples(counts_files)

        # Make the matrix.
        pj = PythonJob('make ' + str(out_matrix_path), [], self.make_counts_matrix,

                       counts_files=counts_files,
                       matrix_out_path=out_matrix_path,
                       reads_to_normalise_to=read_files,
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
                           counts_names: List[str] = None,
                           reads_to_normalise_to: List[Path] = None,
                           add_sites_col: bool = False,
                           normalise_by_site_len: bool = False) -> np.array:
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
                "| awk -v OFS='\\t' {'print $1,$2,$6,$7,$8,$9'}",
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














import logging
import sys
from pathlib import Path
from threading import Thread
from time import sleep
from typing import Dict, Callable, Union, Tuple, List

from components.deeptools import generate_bed_matrix
from utils.fetch_files import get_matching_files
from utils.job_formatter import JobBuilder, PythonJob, ExecParams
from utils.job_manager import JobManager
from utils.mark_colour_map import get_colour
from utils.path_manager import cmdify, PathManager
from utils.utils import TargetedDesign, ExperimentalDesign, vis_namesafe

logger = logging.getLogger(__file__)


def add(dic, key, val):
    if key in dic:
        dic[key].append(val)
    else:
        dic[key] = [val]

def merge(dic, dic2):
    for key in dic2:
        if key in dic:
            dic[key].extend(dic2[key])
        else:
            dic[key] = dic2[key]

NameSafeDict = Dict[Tuple[str, Path], List[str]]

class HeatmapBuilder:
    """

    Programmatically builds heatmaps. Contains numerous safeguards to ensure
    that the final heatmaps are not mislabelled.

    === Private Attributes ===

    job_manager: The job manager to use when executing jobs.
    job_builder: The job manager to use when building jobs.
    file_manager: Used to manage and navigate the project directory.
    start_array: Command to call in order to start a job array.
    stop_array: Command to call in oder to stop a job array
    design: The experimental design to be used.
    """

    _job_manager: JobManager
    _job_builder: JobBuilder
    _file_manager: PathManager
    _start_array: Callable
    _stop_array: Callable

    def __init__(self, design: ExperimentalDesign, job_manager: JobManager,
                 job_builder: JobBuilder, file_manager: PathManager,
                 start_array: Callable, stop_array: Callable) -> None:
        """
        :param job_manager: The job manager to use when executing jobs.
        :param job_builder: The job manager to use when building jobs.
        :param file_manager: Used to manage and navigate the project directory.
        :param start_array: Command to call in order to start a job array.
        :param stop_array: Command to call in oder to stop a job array
        :param design: The experimental design to be used.
        """
        self._design = design
        self._job_manager = job_manager
        self._job_builder = job_builder
        self._file_manager = file_manager
        self._start_array = start_array
        self._stop_array = stop_array

    def gen_heatmap(self,
                    name: str,
                    samples: List[str],
                    bigwigs: Path,
                    peaks: Union[Path, List[Tuple[str, Path]]],
                    scale_factor: Union[int, List[int], str],
                    out_dir: Path,
                    cores_per_job: int,
                    namesafe_check: bool,
                    colour: str,
                    verbose: bool = False,
                    plotheatmap_args: str =
                    '--sortRegions descend --sortUsing mean') -> None:
        """
        Generates a heatmap from provided peak and bigwig files.

        This function creates a heatmap visualization using specified bigwig files
        and peak information. It supports either a single merged peak file or a list
        of annotated peak files. The function also performs validation checks on
        file names and logs any violations. The generated heatmap is saved as an SVG file.

        Parameters:
        name (str): The name of the heatmap.
        samples (List[str]): A list of sample names to be matched with bigwig files.
        bigwigs (Path): Path to the directory containing bigwig files.
        peaks (Union[Path, List[Tuple[str, Path]]]): Either a single Path to a merged
            peak file or a list of tuples containing annotated peak information.
        scale_factor (int): Scale factor for the heatmap intensity.
        out_dir (Path): Output directory path where results will be saved.
        cores_per_job (int): Number of cores to allocate per job.
        namesafe_check (bool): Flag to enable or disable filename validation checks.
        colour (str): Colour for the heatmap.
        verbose (bool, optional): Flag to enable verbose logging. Defaults to False.

        Returns:
        None: This function does not return a value but saves the heatmap image in the output directory.

        Raises:
        ValueError: If the provided `peaks` parameter does not conform to the expected types.

        Example:
        >>> heatmap_generator.gen_heatmap("sample_heatmap", ["sample1", "sample2"],
                                          Path("/path/to/bigwigs"),
                                          Path("/path/to/peaks.bed"), 5,
                                          Path("/output/directory"), 2, True,
                                          "red", True)

        This example will generate a heatmap named 'sample_heatmap' using bigwig files
        corresponding to 'sample1' and 'sample2' in the specified directory, along with
        the peaks information from 'peaks.bed', and save it in '/output/directory'.
        """
        violation_string = ''

        # Fetch Peaks and bigwigs.
        bigwigs = get_matching_files(bigwigs, samples, filetype='bw',
                                     under_delim=True, paths=True,
                                     one_to_one=True)

        if isinstance(self._design, TargetedDesign):
            cond = 'treatment'
        else:
            cond = 'condition'

        treatments = [self._design.query(get=cond,
                                         filters={'sample_name': sn})[0]
                      for sn in samples]

        # Handle creation of Rows
        # Compile all peaks into a single merged peaksfiles
        if isinstance(peaks, Path):
            merged_peaks = out_dir / f'{name}_common_peaks.bed'

            row_names, row_desc, namesafe_vio2 = self.get_merged_peaks(
                peaks, samples, merged_peaks, name, verbose)

            if namesafe_vio2:
                violation_string += 'Peak Violations\n'
                violation_string += vis_namesafe(namesafe_vio2)
            row_names = [row_names]
            peaks = [merged_peaks]

        # Annotate a list of peak files to use as the rows.
        elif isinstance(peaks, list):
            peaks, row_names, row_desc = self.get_annotated_peaks(peaks)

        else:
            raise ValueError('mark_to_peaks contains invalid dictionary '
                             'values.')

        # Perform check to ensure that correct labels and files have been used.
        if namesafe_check:
            new_violations = self._design.namesafe_check(samples, bigwigs)
            if new_violations:
                violation_string += 'Bigwig Violations\n'
                violation_string += vis_namesafe(new_violations)
                logger.info(violation_string)

        if verbose and violation_string:
            print(violation_string, file=sys.stderr)
        if violation_string:
            logger.info(violation_string)

        # Make string representation of the figure.
        col_strs = ''
        for i, bw in enumerate(bigwigs):
            filename = bw.name
            col_name = treatments[i]
            col_strs += '\n\n\t' + col_name + '\t' + '-' + '\t' + filename

        s = f'\t===   Figure Configuration for {name}   ===\n' \
            f'\n' \
            f'Columns:' \
            f'{col_strs}\n' \
            + row_desc

        if verbose:
            print(s)
        logger.info(s)

        # Make matrix.
        matrix = out_dir / f'{name}.gz'
        pj = PythonJob('make ' + str(matrix),
                       [],
                       generate_bed_matrix,

                       beds=peaks,
                       bigwigs=bigwigs,
                       column_names=treatments,
                       row_names=row_names,
                       out_path=matrix,
                       jobs=self._job_manager,
                       builder=self._job_builder,
                       start_array=self._start_array,
                       stop_array=self._stop_array,
                       path_manager=self._file_manager,
                       num_cores=cores_per_job
                       )
        self._job_manager.execute_lazy(pj)

        img = out_dir / f'{name}.svg'

        arrayable = \
            ExecParams(max_runtime=(0, 1, 0), num_cores=cores_per_job,
                       ram_per_core=1024 * 4, builder=self._job_builder,
                       wait=False)
        if isinstance(scale_factor, list):
            scale_factor = ','.join(scale_factor)
        cmd = cmdify(
            "plotHeatmap",
            "-m", matrix,
            "-o", img,
            "--yMin 0 --zMin 0",
            "--yMax", scale_factor,
            "--zMax", scale_factor,
            "--colorList", colour, plotheatmap_args
        )
        self._job_manager.execute_lazy(cmd, arrayable)

    def get_annotated_peaks(self, peaks: List[Tuple[str, Path]]) -> \
            Tuple[List[Path], List[str], str]:
        """
        :param peaks: A list of peaks and the name of the row to which they
            corespond.
        :return: Paths to each of the peak files, the annotated column names,
            a repesentation of the rows added.
        """

        bedfiles = []
        row_names = []
        row_desc = f'Rows:\n\n'

        for rowname, bedfile in peaks:
            bedfiles.append(bedfile)

            with open(bedfile) as peaksfile:
                n_merged_peaks = len(peaksfile.readlines())
            row_name = f'{rowname} (n = {n_merged_peaks})'
            row_names.append(row_name)
            row_desc += f'\t{row_name}\t-\t{bedfile.name}\n\n'

        return bedfiles, row_names, row_desc

    def get_merged_peaks(self, directory: Path, sample_names: List[str],
                         out_peakfile: Path, merged_bedname: str,
                         verbose: bool = False) \
            -> Tuple[str, str, NameSafeDict]:
        """

        :param directory: Directory containing the peaks to merge.
        :param mark: Beds associated with the given mark will be merged.
        :param verbose: Display additonal information.
        :param out_peakfile: Where to place the merged peaks file.
        :return: The row name, the description of the row name, and any
            namesafe violations found during the merge process.
        """
        peaks = get_matching_files(directory, sample_names,
                                   filetype='(((narrow|broad)Peak)|peaks.bed)',
                                   under_delim=True, paths=True,
                                   one_to_one=True)

        self._job_manager.execute_lazy(
            cmdify('cat', *peaks,
                   '| bedtools sort',
                   '| bedtools merge',
                   '>', out_peakfile))

        with open(out_peakfile) as peaksfile:
            n_merged_peaks = len(peaksfile.readlines())
        row_name = f'All Merged {merged_bedname} Peaks (n = {n_merged_peaks})'

        namesafe_violations = self._design.namesafe_check(sample_names, peaks)

        # Display merging process.
        merge_str = ''
        merge_str += '\t --- Peak Merging ---\n'
        for file in peaks:
            with open(file, 'r') as peakfile:
                sub_n_peaks = len(peakfile.readlines())
            merge_str += '\t' + str(sub_n_peaks) + '\t' + file.name + '\n'
        merge_str += f'Final Merged:\n\t{n_merged_peaks}\t' \
                         f'{out_peakfile.name}\n'
        if verbose:
            print(merge_str)
        logger.info(merge_str)

        row_desc = f'Rows:\n\n' \
                   f'\t{row_name}\t-\t{out_peakfile.name}\n\n'

        # Return path to merged peaks, name of merged peaks, and a description
        # to be used as a row name.
        return row_name, row_desc, namesafe_violations

    def heatmaps_by_mark(self,
                         mark_to_bigwigs_dir: Dict[str, Path],
                         mark_to_peaks: Dict[str, Union[Path,
                                                  List[Tuple[str, Path]]]],
                         mark_to_scale_factor: Dict[str, int],
                         out_dir: Path,
                         cores_per_job: int = 1,
                         verbose: bool = False,
                         namesafe_check: bool = True,
                         run_async: bool = False,
                         plotheatmap_args: str =
                         '--sortRegions descend --sortUsing mean'
                         ) -> None:
        """
        Generate heatmaps for each mark.

        :param cores_per_job:
        :param mark_to_bigwigs_dir: Maps a mark to the directory containing the
            bigwigs to use for that mark.
        :param mark_to_peaks:
            Maps a mark to a directory containing bed files or a list of
            bedfiles alongside their row names.

            Iff the dict maps to a path, extract bedfiles for each treatment for
            the given mark and merge them into a single common peaks bedfile.

            Iff the dict maps to a list of Tuple[str, Path], use each of the Paths
            as the row bedfile and the name as the name of that column.

        :param mark_to_scale_factor: Provides a scaling factor for each of the
            given marks.
        :param out_dir: The directory into which all the output figures and
            intermediary files will be placed.
        :param verbose: Print addiitonal information regarding the process to the
            console.
        :param namesafe_check: Perform addition checks to ensure that the correct
            files are being used independent of specifications provided by
            <design>.
        :return: None.
        """
        jobs = []

        # Generate matrices.
        for mark in set(self._design.query('target')):
            # Configure output directory
            mark_common_dir = self._file_manager.make(out_dir / mark)
            samples = self._design.query('sample_name', {'target': mark})

            jobs.append(PythonJob(
                f'Make {mark_common_dir}', [],

                to_execute=self.gen_heatmap,
                name=mark,
                samples=samples,
                bigwigs=mark_to_bigwigs_dir[mark],
                peaks=mark_to_peaks[mark],
                scale_factor=mark_to_scale_factor[mark],
                out_dir=mark_common_dir,
                cores_per_job=cores_per_job,
                namesafe_check=namesafe_check,
                colour='white,' + get_colour(mark),
                verbose=verbose,
                plotheatmap_args=plotheatmap_args

            ))

        threads = []
        for job in jobs:
            threads.append(Thread(target=self._job_manager.execute,
                                  args=[job]))
        for thread in threads:
            thread.start()

        if not run_async:
            for thread in threads:
                thread.join()






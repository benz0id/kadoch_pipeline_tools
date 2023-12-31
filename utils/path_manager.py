import os
import shutil
import sys
from datetime import datetime
from pathlib import Path
import logging
from typing import List, Union
from constants.data_paths import PIPELINE_BACKUP_DIR

EXECUTABLES_TO_SAVE = [
            'Rmd',
            'R',
            'sh',
            'py'
        ]

global wd

def quotes(s: Union[str, Path]) -> str:
    """
    Wraps the given string in quotes. Useful when placing filepaths into commands.
    :param s: The string to be wrapped in quotes.
    :return: <s>, wrapped in quotes.
    """
    if isinstance(s, Path):
        s = str(s)

    return f'"{s}"'


def cmdify(*args: Union[str, Path, int, float]) -> str:
    """
    Attempts to format the given args as strings and join them into a command.
    :param args: Some number of objects.
    :return:
    """
    global wd
    formatted = []
    for arg in args:
        if isinstance(arg, str):
            formatted.append(arg)
        elif isinstance(arg, Path):
            if wd in arg.parents:
                arg = os.path.relpath(arg, wd)
            formatted.append(quotes(arg))
        elif isinstance(arg, int):
            formatted.append(str(arg))
        elif isinstance(arg, float):
            formatted.append(str(arg))
        else:
            raise ValueError(f'Unrecognised type: {repr(arg)}')

    return ' '.join(formatted)


def sizeof_fmt(num: float) -> str:
    """
    from https://stackoverflow.com/questions/1094841/get-human-readable-version-of-file-size

    Get a human-readable representation of the given amount of data in bytes.

    :param num: The amount of data, in bytes.
    :return: A humanreadable representation of data size.
    """
    for unit in ("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"):
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}B"
        num /= 1024.0
    return f"{num:.1f}YiB"


def get_all_outside_dir(directory: Path, files: List[Path]) -> List[Path]:
    """
    Returns all <files> that do not lie within <directory>.
    :param directory: The containing directory.
    :param files: Files to be checked for presence in directory.
    :return: A list of all <files> not in <directory>
    """
    not_in = []
    for file in files:
        if directory not in file.parents:
            not_in.append(file)
    return not_in


def get_total_size(files: List[Path]) -> int:
    """
    Returns the total size of all the given files in bytes.
    :param files: A list of filepaths leading to real files.
    :return: The size of all <files> in bytes.
    """
    agg_size = 0
    for file in files:
        agg_size += os.path.getsize(file)
    return agg_size


class PathManager:
    """

    === Description ===
    Responsible for managing the structure of essential project components.

    === Basic Project Structure ===

    project_dir             (self.project_dir)
    |---- pipeline.py
    |---- sample_sheets     (self.sample_sheets_dir)
    |---- sequencing_dir    (self.sequencing_dir)
    |     |---- <your sequencing results>
    |     |---- demult_stats(self.demult_stats)
    |     |---- fastqs      (self.fastqs_dir)
    |     |---- fastqc      (self.fastqc_dir)
    |---- cache_dir (self.cache_dir)
    |     |---- <any caches you create with cache manager>
    |---- logs              (self.logging_directory)
    |     |---- info_archive
    |     |---- debug_archive
    |     |---- general_archive
    |     |---- INFO.log
    |     |---- DEBUG.log
    |---- .purgeable        (self.purgeable_files_dir)
    |---- .slurm_scripts
    |---- sbatch.out
    |---- sbatch.err

    === Brief Descriptions ===

    project_dir
        The main project directory.
    pipeline.py
        The pipeline that you have written.
    sample_sheets
        Directory containing sample sheets.
    sequencing_dir
        Directory containing raw sequencing data.
    <your sequencing results>
        Raw Illumina sequencing data.
    demult_stats
        Results produced by demultiplexing.
    fastqs
        Raw fastqs generated after running basecalling.
    fastqc
        Results from fastqc analysis.
    cache_dir
        Directory containing text files that log previously executed commands
        managed by the <cache_manager>.
    logs
        Contains various logs.
    .info_archive
        All previously generated info logs.
    .debug_archive
        All previously generated debug logs.
    .general_archive
        Misc logs. Also contains backups of previous versions of pipeline.py.
    INFO.log
        Info log from last execution.
    DEBUG.log
        Debug log from last execution.
    .purgeable
        Temporary files that can be removed after project completion.
    .slurm_scripts
        Directory for storing temporary sbatch scripts. Managed by Slurmifier
    sbatch.out
        Log file for all sbatch output.
    sbatch.err
        Log file for all sbatch err.
    """
    project_dir: Path
    cache_dir: Path
    logging_directory: Path
    info_archive: Path
    debug_archive: Path
    general_archive: Path
    home_dir: Path
    sequencing_dir: Path
    purgeable_files_dir: Path
    sample_sheets_dir: Path
    fastqs_dir: Path
    sample_sheet_path: Path
    sequencing_results: Path
    fastqc_dir: Path
    demult_stats: Path

    def __init__(self, working_dir: Path, verbose: bool = False) -> None:
        self.verbose = verbose
        self.project_dir = working_dir
        global wd
        wd = working_dir
        self.configure_required_dirs()
        self.backup_scripts()

    def backup_scripts(self) -> None:
        """
        Copies a backup of all scripts in the working directory to the backup
        directory.
        :return:
        """

        time_str = datetime.now().strftime("%m-%d-%Y_%I:%M:%S_%p")
        files = os.listdir(self.project_dir)
        for file in files:
            filepath = self.project_dir / file
            filetype = '.'.join(file.split('.')[1:])
            raw_filename = file.split('.')[0]

            p1 = filepath.parent
            p2 = p1.parent

            out_dir = self.make(PIPELINE_BACKUP_DIR / p2.name)
            out_dir = self.make(out_dir / p1.name)
            out_name = raw_filename + time_str + '.' + filetype

            if filetype in EXECUTABLES_TO_SAVE:
                shutil.copy(filepath, out_dir / out_name)

    def add_basic_info(self, sequencing_results_name: str,
                       sample_sheet_name: str) -> None:
        """
        Assuming that the <sequencing_results_name> is a subdirectory of
        self.sequencing_dir and <sample_sheet_name> is within
        self.sample_sheets_dir, stores paths to these folders.

        :param sequencing_results_name: Name of the raw sequencing results dir,
        eg: 230915_NB551325_0672_AHHK2NBGXT_GX11457.

        :param sample_sheet_name: The name of the raw sample sheet.

        === Example ===

        >>> pm = PathManager(Path('/project_dir'))
        >>> pm.add_basic_info("230915_NB551325_0672_AHHK2NBGXT_GX11457",
        ...                   "My_Sample_Sheet.csv")
        """
        self.sample_sheet_path = self.sample_sheets_dir / sample_sheet_name
        self.sequencing_results = self.sequencing_dir / sequencing_results_name

    @staticmethod
    def safe_make(dir: Path) -> bool:
        """Makes the given directory if it does not already exist.

        :return: Whether a directory was created."""
        exists = dir.exists()
        if not exists:
            os.makedirs(dir)
        return not exists

    def rel_make(self, relative_path: Union[Path, str]) -> Path:
        """Creates the given directory in the working directory."""
        if isinstance(relative_path, str):
            relative_path = Path(relative_path)
        full_path = self.project_dir / relative_path
        self.safe_make(full_path)
        return full_path

    def make(self, new_directory_path: Path) -> Path:
        """
        Makes a new directory iff it does not already exist, returns the given
         path.
        :param new_directory_path: Path to the directory to be created.
        :return: Path to the directory to be created.
        """
        self.safe_make(new_directory_path)
        return new_directory_path

    def configure_required_dirs(self) -> None:
        """Configure the directories required for pipeline utility functionality."""
        self.cache_dir = self.project_dir / '.pipeline_cache'
        self.safe_make(self.cache_dir)

        self.logging_directory = self.project_dir / 'logs'
        self.safe_make(self.logging_directory)

        self.purgeable_files_dir = self.project_dir / '.purgeable'
        self.safe_make(self.purgeable_files_dir)

        self.home_dir = Path(os.path.expanduser(""))

        self.info_archive = self.logging_directory / '.info_archive'
        self.safe_make(self.info_archive)

        self.debug_archive = self.logging_directory / '.debug_archive'
        self.safe_make(self.debug_archive)

        self.general_archive = self.logging_directory / '.general_archive'
        self.safe_make(self.general_archive)

        self.sequencing_dir = self.project_dir / 'sequencing'
        self.safe_make(self.sequencing_dir)

        self.fastqs_dir = self.sequencing_dir / 'fastqs'
        self.safe_make(self.fastqs_dir)

        self.fastqc_dir = self.sequencing_dir / 'fastqc'
        self.safe_make(self.fastqc_dir)

        self.sample_sheets_dir = self.project_dir / 'sample_sheets'
        self.safe_make(self.sample_sheets_dir)

        self.demult_stats = self.sequencing_dir / 'demult_stats'
        self.safe_make(self.demult_stats)

        self.move_logs_to_archive()
        self.local_pipeline_backup()
        self.configure_logging()

    def local_pipeline_backup(self) -> None:
        """
        Make a backup of all python files in project directory.
        :return:
        """
        
        time_str = datetime.now().strftime("%m-%d-%Y_%I:%M:%S_%p")
        for filename in os.listdir(self.project_dir):
            filetype = str(filename).split('.')[-1]
            for ident in EXECUTABLES_TO_SAVE:
                if filetype == ident:
                    os.system(cmdify(
                        'cp', self.project_dir / filename,
                              self.general_archive /
                              (filename[-3] + time_str + '.' + filetype)))

    def move_logs_to_archive(self) -> None:
        """
        Move any logfiles in logging directory to archive.
        """
        existing_logfiles = []
        for filepath in os.listdir(self.logging_directory):
            if '.log' in str(filepath):
                existing_logfiles.append(filepath)

        def targeted_archive(identifier: str, archive_file: Path):
            to_rem = []
            for i, logfile in enumerate(existing_logfiles):
                if identifier in str(logfile):
                    shutil.move(self.logging_directory / logfile,
                                archive_file / os.path.basename(logfile))
                    to_rem.append(i)

            for i in sorted(to_rem, reverse=True):
                existing_logfiles.pop(i)

        targeted_archive('DEBUG', self.debug_archive)
        targeted_archive('INFO', self.info_archive)
        targeted_archive('', self.general_archive)

    def configure_logging(self) -> None:
        """Configures a simple logger."""
        format = '%(asctime)s - %(levelname)s: %(message)s'
        datefmt = '%m/%d/%Y %I:%M:%S %p'

        formatter = logging.Formatter(fmt=format, datefmt=datefmt)

        info_file = self.logging_directory / datetime.now().strftime(
            "INFO-%m-%d-%Y_%I:%M:%S_%p.log")
        debug_file = self.logging_directory / datetime.now().strftime(
            "DEBUG-%m-%d-%Y_%I:%M:%S_%p.log")

        logging.root.handlers = []

        info_handler = logging.FileHandler(info_file)
        debug_handler = logging.FileHandler(debug_file)

        info_handler.setLevel(logging.INFO)
        debug_handler.setLevel(logging.DEBUG)
        info_handler.setFormatter(formatter)
        debug_handler.setFormatter(formatter)

        logging.root.addHandler(info_handler)
        logging.root.addHandler(debug_handler)
        logging.root.setLevel(0)
        if self.verbose:
            logging.root.addHandler(logging.StreamHandler(sys.stdout))

import os
import shutil
import sys
from datetime import datetime
from pathlib import Path
import logging
from typing import List, Union


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
    formatted = []
    for arg in args:
        if isinstance(arg, str):
            formatted.append(arg)
        elif isinstance(arg, Path):
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
    Class responsible for the creation and management of paths within the project directory.
    Also, currently handles logging.

    === Public Attributes ===

    working_dir: The working directory of the project.

    """
    project_dir: Path
    cache_filepath: Path
    logging_directory: Path
    info_archive: Path
    debug_archive: Path
    general_archive: Path
    home_dir: Path

    def __init__(self, working_dir: Path, verbose:bool = False) -> None:
        self.verbose = verbose
        self.project_dir = working_dir
        self.configure_required_dirs()

    @staticmethod
    def safe_make(dir: Path) -> bool:
        """Makes the given directory if it does not already exist.

        :return: Whether a directory was created."""
        exists = dir.exists()
        if not exists:
            os.makedirs(dir)
        return not exists

    def make_dir(self, relative_path: Path) -> Path:
        """Creates the given directory in the working directory."""
        full_path = self.project_dir / relative_path
        self.safe_make(full_path)
        return full_path

    def configure_required_dirs(self) -> None:
        """Configure the directories required for pipeline utility functionality."""
        self.cache_filepath = self.project_dir / '.pipeline_cache'
        self.safe_make(self.cache_filepath)

        self.logging_directory = self.project_dir / 'logs'
        self.safe_make(self.logging_directory)

        self.home_dir = Path(os.path.expanduser(""))

        self.info_archive = self.logging_directory / '.info_archive'
        self.safe_make(self.info_archive)

        self.debug_archive = self.logging_directory / '.debug_archive'
        self.safe_make(self.debug_archive)

        self.general_archive = self.logging_directory / '.general_archive'
        self.safe_make(self.general_archive)

        self.move_logs_to_archive()
        self.pipeline_backup()
        self.configure_logging()

    def pipeline_backup(self) -> None:
        """
        Make a backup of all python files in project directory.
        :return:
        """
        time_str = datetime.now().strftime("%m-%d-%Y_%I:%M:%S_%p.py")
        for filename in os.listdir(self.project_dir):
            if '.py' in filename:
                os.system(cmdify(
                    'cp', self.project_dir / filename,
                          self.general_archive / (filename[-3] + time_str)))

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

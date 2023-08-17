import os
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


def config_root_logger(logging_dir: Path) -> None:
    """Configures a simple logger."""
    filename = datetime.now().strftime("%m/%d/%Y_%I:%M:%S_%p.log")
    logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
                        filename=logging_dir / filename, level=0)
    logging.root.handlers[0].setLevel(0)


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
        if file not in directory:
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

    === Public Attributes ===

    working_dir: The working directory of the project.

    """
    project_dir: Path
    cache_filepath: Path
    logging_directory: Path

    def __int__(self, working_dir: Path, ) -> None:
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
        config_root_logger(self.logging_directory)

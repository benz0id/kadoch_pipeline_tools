import logging
import os
import shutil
import sys
from pathlib import Path
from typing import List, Union

from utils.job_formatter import Job
from utils.path_manager import PathManager, get_all_outside_dir, get_total_size, sizeof_fmt

logger = logging.getLogger(__name__)


class CacheManager:
    """

    === Description ===
    Class designed to reduce the needless repetition of compute-intensive tasks and allow for removal of redundant
    or otherwise unneeded data.

    === Private Attributes ===

    path_manager: PathManager to guide placement of files.

    cur_commands: A list of all commands that have been in the current execution of the program.

    past_cmds: A list of commands run during previous executions of the program.

    purgeable_data: Data to be purged after successful pipeline execution.

    strict: Changes class behaviour such that if a change in a pipeline step is detected, all downstream steps are
        removed from the cache and repeated. Use to ensure that no downstream steps are run using outdated data.

        ***WARNING*** Any modification or addition of ANY upstream step that is cached by this CacheManager will result
        in the re-execution of all subsequent cached commands if this parameter is True.

    name: The name of this cache manager. Allows for multiple independent caches to function in a single pipeline.
    """
    _name: str
    _path_manager: PathManager

    _cur_exec_cmds: List[str]
    _cur_skipped_cmds: List[str]
    _past_cmds: List[str]

    _cache_path: Path

    _strict: bool
    _pipeline_integrity_maintained: bool

    _purgeable_data: List[Path]

    def __init__(self, path_manager: PathManager, name: str = 'basic_cache', strict: bool = False) -> None:
        self._path_manager = path_manager

        self._strict = strict
        self._pipline_integrity_maintained = True

        self._cache_path = path_manager.cache_filepath / name
        self._cache_path.touch()
        self._read_cache()
        self._cur_exec_cmds = []
        self._cur_skipped_cmds = []
        self._purgeable_data = []
        self._name = name

    def _read_cache(self) -> None:
        """Read previously executed commands from the cache."""
        with open(self._cache_path, 'r') as file:
            lines = file.readlines()
            cmds = ''.join(lines)
            cmds = cmds.split('\n<END OF COMMAND>\n')
            self._past_cmds = cmds

        logger.debug('Previously logged commands:\n\t\t' + '\n\t\t'.join(self._past_cmds))

    def do_execution(self, job: Job) -> bool:
        """
        Returns whether this job has already been run in this pipeline.

        :param Job job: A job that may or may not have been run in previous iterations of this pipeline.
        """
        return job.get_cmd() not in self._past_cmds

    def cache_skipped(self, job: Job) -> None:
        """
        Store the fact that a job was skipped.
        :param job:
        :return:
        """
        self._cur_skipped_cmds.append(job.get_cmd())

    def write_cache(self) -> None:
        """Writes the currently stored cache to a file."""
        with open(self._cache_path, 'w') as f:
            f.writelines([cmd + '\n<END OF COMMAND>\n' for cmd in self._past_cmds] + ['\n'])

    def cache_execution(self, job: Job) -> None:
        """
        Updates cache with job information if necessary.

        :param job: The job that executed the command to be stored.
        """
        self._cur_exec_cmds.append(job.get_cmd())

        if job.get_cmd() not in self._past_cmds:
            logger.info(f'Detected new step in cached pipeline: {self._name}\n\t\t New step: {job.get_cmd()}.')
            self._past_cmds.append(job.get_cmd())

            # Pipeline modification detected!
            if self._strict and self._pipline_integrity_maintained:
                self._pipline_integrity_maintained = False
                commands_to_be_repeated = list(set(self._past_cmds) - set(self._cur_exec_cmds))

                logger.info(f"{self._name} is configured to be strict. As such, the following steps will be repeated"
                            f" if they are still in the pipeline:\n\t\t" + '\n\t\t'.join(commands_to_be_repeated))
                self._past_cmds = self._cur_skipped_cmds[:] + self._cur_exec_cmds[:]

            self.write_cache()

    def add_purgeable_data(self, purgeable: Union[List[Path], Path]) -> None:
        """
        Adds all the given filepaths to the list of files that can be purged without consequence after the
        execution of the pipeline.

        :param List[Path] purgeable: A list of files that will be deleted iff self.purge_data is called.
        """
        if isinstance(purgeable, list):
            self._purgeable_data.extend(purgeable)
            logger.debug(
                'Adding the following files to be purged.\n\t\t' + '\n\t\t'.join([str(path) for path in purgeable]))
        else:
            self._purgeable_data.append(purgeable)
            logger.debug(
                'Adding the following file to be purged.\n\t\t' + str(purgeable))

    def purge_data(self, req_user_input: bool = True, size_limit: int = 20) -> None:
        """
        Deletes all files in self._purgeable_files. All files must be contained within the current directory.

        :param req_user_input: Whether to prompt the user before deleting the given files.
        :param size_limit: The absolute size limit of all files to be deleted, in GiB.
        """
        if not self._purgeable_data:
            print('No data to purge.')
            logger.info('No data to purge.')
            return

        to_be_deleted_str = '\n\t\t'.join([str(p) for p in self._purgeable_data])

        logger.info('=== Begging File Purge Process ===')
        logger.info('Files to be purged:\n\t\t' + to_be_deleted_str)

        # Check that all files are in project directory.
        files_not_in_dir = get_all_outside_dir(self._path_manager.project_dir, self._purgeable_data)
        if len(files_not_in_dir) == 0:
            logger.info('All files to be purged are in project directory.')
        else:
            s = 'Purge cannot be completed, as there are several files to be purged that lie outside of the project' \
                'directory, including: \n\t\t' + '\n\t\t'.join([str(s) for s in files_not_in_dir])
            print(s, file=sys.stderr)
            logger.error(s)
            return

        agg_size = get_total_size(self._purgeable_data)
        max_size = size_limit * 1024 ** 3
        if agg_size <= max_size:
            logger.info(f'File size found to be within purge thresholds: '
                        f'{sizeof_fmt(agg_size)} <= {sizeof_fmt(max_size)}')
        else:
            s = f'Purge cannot be completed, as aggregate size of files to be deleted is larger that allowable ' \
                f'amount: {sizeof_fmt(agg_size)} > {sizeof_fmt(max_size)}'
            print(s, file=sys.stderr)
            logger.error(s)
            return

        if req_user_input:
            logger.info('Requesting user authorization for file purge.')

            print("=== Manual Approval for File Purge ===")
            print("Files to be deleted: \n\t\t" + to_be_deleted_str)
            print("Total size: ", sizeof_fmt(agg_size))

            prompt = 'Approve? (Y/n)'
            usr_in = input(prompt).strip()
            while usr_in != 'Y' and usr_in != 'n':
                print("Invalid input. Try again...")
                usr_in = input(prompt).strip()

            if usr_in == 'Y':
                s = 'Purge approved. Beginning purge.'
                logger.info(s)
                print(s)
            else:
                s = "Purge rejected by user."
                logger.info(s)
                print(s)
                return

        # Finally, remove files.
        while self._purgeable_data:
            filepath = self._purgeable_data.pop()
            if filepath.is_dir():
                shutil.rmtree(filepath)
            else:
                os.remove(filepath)












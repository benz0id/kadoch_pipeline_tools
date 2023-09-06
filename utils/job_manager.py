import logging
from pathlib import Path
from typing import List, Union
from utils.cache_manager import CacheManager
from utils.path_manager import PathManager
from utils.job_formatter import ExecParams, Job, JobBuilder

logger = logging.getLogger(__name__)

default_job = ExecParams(max_runtime=(0, 0, 10), num_cores=1,
                             ram_per_core=128, builder=JobBuilder())
class JobManager:
    """
    Responsible for managing the execution of jobs.

    === Private Attributes ===
    cache_manager: CacheManager responsible for directing this JobManagers execution decisions.
    path_manager: PathManager to guide placement of files.
    prevent_execution: Whether to actually execute commands.
    """

    def __init__(self, cache_manager: CacheManager, path_manager: PathManager,
                 verbose: bool = True):
        self._cache_manager = cache_manager
        self._path_manager = path_manager
        self._prevent_execution = False
        self._verbose = verbose

    def disable_execution(self) -> None:
        """
        Configures this JobManger as a dummy manager that des not actually execute commands.
        """
        self._prevent_execution = True

    def enable_execution(self) -> None:
        """
        Configures this JobManger as a real manager that actually executes commands.
        :return:
        """
        self._prevent_execution = False

    def execute(self, cmd: Union[str, Job],
                exec_params: ExecParams = default_job) -> None:
        """Execute the stored in <cmd>, outputting to stderr when possible.

        :param str cmd: The command/job to be executed.
        :param ExecParams exec_params: Parameters used to execute the command.
        """
        if isinstance(cmd, str):
            job = exec_params.builder.prepare_job(cmd, exec_params)
        else:
            job = cmd

        if self._verbose and not self._prevent_execution:
            print(repr(exec_params), '\n\t/bin/bash', job.get_cmd())
        if not self._prevent_execution:
            job.execute()

    def execute_lazy(self, cmd: Union[str, Job],
                     exec_params: ExecParams = default_job) -> bool:
        """Executes <cmd> iff it has not been called this run and all previous steps in the pipeline have been
        successfully executed lazily.

        :param exec_params: The parameters used to configure the command.
        :param str cmd: The command to be executed.
        :return: Whether cmd was executed.
        """
        if isinstance(cmd, str):
            job = exec_params.builder.prepare_job(cmd, exec_params)
        else:
            job = cmd

        execute_job = self._cache_manager.do_execution(job)
        if execute_job:
            if not self._prevent_execution:
                logger.info(f'Lazily executing `{cmd}`.')
                job.execute()
            self._cache_manager.cache_execution(job)
            logger.info('Job complete. Storing command in cache.')
            return True
        else:
            logger.info(f'Skipping `{cmd}`, as it already exists in the cache.')
            self._cache_manager.cache_skipped(job)
            return False

    def execute_purgeable(self, cmd: Union[str, Job],
                          purgeable: Union[Path, List[Path]],
                          exec_params: ExecParams = default_job,
                          lazy: bool = True) -> bool:
        """Executes <cmd> iff it has not been called this run and all previous steps in the pipeline have been
        successfully executed lazily. Adds the given files to the list of files that can be purged after the pipeline
        has sucessfully finished running.

        :param exec_params: The parameters used to configure the command.
        :param lazy: Whether to execute the command lazily.
        :param str cmd: The command to be executed.
        :param purgeable: A list of file paths to be deleted after the pipeline has completed.
        """
        logger.info('Command received with output to be purged.')

        if lazy:
            executed = self.execute_lazy(cmd, exec_params)
        else:
            executed = True
            self.execute(cmd, exec_params)

        if executed:
            self._cache_manager.add_purgeable_data(purgeable)

        return executed



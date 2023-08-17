from pathlib import Path
from typing import List
from utils.cache_manager import CacheManager
from utils.path_manager import PathManager
from utils.job_formatter import ExecParams


class JobManager:
    """
    Responsible for managing the execution of jobs.

    === Private Attributes ===
    cache_manager: CacheManager responsible for directing this JobManagers execution decisions.
    path_manager: PathManager to guide placement of files.
    """

    def __init__(self, cache_manager: CacheManager, path_manager: PathManager):
        self.cache_manager = cache_manager
        self.path_manager = path_manager

    def execute(self, cmd: str, exec_params: ExecParams) -> None:
        """Execute the stored in <cmd>, outputting to stderr when possible.

        :param str cmd: The command to be executed.
        :param ExecParams exec_params: Parameters used to execute the command.
        """

        job = exec_params.builder.prepare_job(cmd, exec_params)
        job.execute()

    def execute_lazy(self, cmd: str, exec_params: ExecParams) -> None:
        """Executes <cmd> iff it has not been called this run and all previous steps in the pipeline have been
        successfully executed lazily.

        :param str cmd: The command to be executed.
        """
        job = exec_params.builder.prepare_job(cmd, exec_params)
        execute_job = self.cache_manager.do_execution(job)
        if execute_job:
            job.execute()
            self.cache_manager.cache_execution(job)

    def execute_purgeable(self, cmd: str, purgeable: List[Path], exec_params: ExecParams) -> None:
        """Executes <cmd> iff it has not been called this run and all previous steps in the pipeline have been
        successfully executed lazily. Adds the given files to the list of files that can be purged after the pipeline
        has sucessfully finished running.

        :param str cmd: The command to be executed.
        :param purgeable: A list of file paths to be deleted after the pipeline has completed.
        """
        self.cache_manager.add_purgeable_data(purgeable)

        job = exec_params.builder.prepare_job(cmd, exec_params)
        execute_job = self.cache_manager.do_execution(job)
        if execute_job:
            job.execute()
            self.cache_manager.cache_execution(job)



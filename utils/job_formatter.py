import abc
import collections
import logging
import os
from dataclasses import dataclass
from subprocess import run
from typing import Union, Any, List

from patterns.observer import Observable, Observer
from utils.defaults import MAX_RUNTIME, NUM_CORES, RAM_PER_CORE

Runtime = collections.namedtuple('Runtime', ['days', 'hours', 'minutes'])
logger = logging.getLogger('root')


class Job(Observable):
    """
    === Description ===
    A simple job that stores its execution information.

    === Private Attributes ===
    run_str: The string that this job runs.
    complete: whether this job has been run.
    """
    _cmd: str

    def __init__(self, cmd: str, observers: List[Observer] = None) -> None:
        super().__init__(observers)
        self._cmd = cmd
        self.complete = False

    def get_cmd(self) -> str:
        return self._cmd

    def execute(self) -> None:
        logger.info('Executing ' + self._cmd)
        os.system(self._cmd)
        self.complete = True
        self.notify_observers()


class JobBuilder:
    """
    === Description ===
    A simple job builder that merely creates jobs to be run directly from the command line.
    """

    def prepare_job(self, cmd: str, exec_params: Any) -> Job:
        return Job(cmd)


@dataclass
class ExecParams:
    """
    === Description ===
    Parameters required for specifying how a job should be run.

    num_cores: The number of cores to use for the job.
    max_runtime: The maximum runtime of the job.
    ram_per_core: The memory to reserve per CPU, in MiB.
    builder: builder class required to format the job for execution on the machine.
    """

    max_runtime: Runtime
    num_cores: int
    ram_per_core: int
    builder: JobBuilder

    def __init__(self, max_runtime: Runtime = MAX_RUNTIME,
                 num_cores: int = NUM_CORES,
                 ram_per_core: int = RAM_PER_CORE,
                 builder: Union[JobBuilder, None] = None) -> None:
        self.max_runtime = Runtime(*max_runtime)
        self.num_cores = num_cores
        self.ram_per_core = ram_per_core
        self.builder = builder

import collections
import logging
import os
import sys
from copy import copy
from dataclasses import dataclass
from threading import Thread
from typing import Union, Any, List, Dict

from patterns.observer import Observable, Observer
from constants.defaults import MAX_RUNTIME, NUM_CORES, RAM_PER_CORE

Runtime = collections.namedtuple('Runtime', ['days', 'hours', 'minutes'])
logger = logging.getLogger(__name__)


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


class PythonJob(Job):
    """
    === Description ===
    A job that allows for the execution of a custom command with given arguments.

    === Private Attributes ===
    to_execute: The command to be executed with args and kwargs.
    args: To be passed to <to_execute>.
    kwargs: To be passed to <to_execute<.

    """
    _run_async: bool

    def __init__(self, cmd: str, observers: List[Observer], to_execute, *args,
                 **kwargs) -> None:
        """
        Executes <to_execute> using <args> and <kwargs>.
        :param cmd: A string representation of the command to be executed.
        :param observers: The observers to notify when this job is executed.
        :param to_execute: The function or method descriptor to pass <args> and <kwargs> upon execution.
        :param run_async: Whether to spawn a thread to execute the job.
        :param args: To be passed to <to_execute>.
        :param kwargs: To be passed to <to_execute>.
        """
        super().__init__(cmd, observers)

        if 'run_async' in kwargs and kwargs['run_async']:
            self._run_async = True
            self._thread = Thread(target=to_execute, args=args, kwargs=kwargs)
            del kwargs['run_async']
        else:
            self._run_async = False

        self._to_execute = to_execute
        self.args = args
        self.kwargs = kwargs

    def get_thread(self) -> Thread:
        return self._thread

    def execute(self) -> None:
        logger.info(
            f'Executing {repr(self._to_execute)} with args: {str(self.args)} and kwargs: {str(self.kwargs)}.')
        if self._run_async:
            self._thread.start()
            return
        self._to_execute(*self.args, **self.kwargs)
        self.notify_observers()


class JobBuilder:
    """
    === Description ===
    A simple job builder that creates jobs to be run directly from the
    command line.
    """

    def prepare_job(self, cmd: str, exec_params: Any) -> Job:
        # Add lines to load the required modules.
        if len(exec_params.get_requirements()) != 0:
            load_line = 'module load '
            for requirement in exec_params.get_requirements():
                load_line += f'{requirement}/{exec_params.get_requirements()[requirement]} '
            load_line += '\n'
        else:
            load_line = ''

        return Job(load_line + cmd)


@dataclass
class ExecParams:
    """
    === Description ===
    Parameters required for specifying how a job should be run.

    num_cores: The number of cores to use for the job.
    max_runtime: The maximum runtime of the job.
    ram_per_core: The memory to reserve per CPU, in MiB.
    builder: Builder class required to format the job for execution on the machine.
        Job is just run in active shell session iff None.
    requires: A list of string names of the programs required to run this job.
    """

    max_runtime: Runtime
    num_cores: int
    ram_per_core: int
    builder: JobBuilder
    _requires: Dict[str, str]
    wait: bool

    def __init__(self, max_runtime: Runtime = MAX_RUNTIME,
                 num_cores: int = NUM_CORES,
                 ram_per_core: int = RAM_PER_CORE,
                 builder: Union[JobBuilder, None] = None,
                 requires: Dict[str, str] = None,
                 wait: bool = True) -> None:
        self.wait = wait
        self.max_runtime = Runtime(*max_runtime)
        self.num_cores = num_cores
        self.ram_per_core = ram_per_core
        self.builder = builder
        if requires is None:
            self._requires = {}
        else:
            self._requires = requires

    def __repr__(self):
        s = f'ExecParams(runtime={self.max_runtime}, ' \
            f'cores={self.num_cores}, ' \
            f'ram={self.ram_per_core}, ' \
            f'builder={self.builder}' \
            f'wait={self.wait}'
        if self._requires:
            s += f', requires={self._requires} '
        return s

    def __str__(self):
        s = f'Executed with ' \
            f'\n\truntime={self.max_runtime} ' \
            f'\n\tcores={self.num_cores} ' \
            f'\n\tram={self.ram_per_core} ' \
            f'\n\tbuilder={self.builder} '\
            f'\n\t wait={self.wait}'
        if self._requires:
            s += f' \n\trequires={self._requires}'
        return s

    def add_requirements(self, new_requirements: Dict[str, str]) -> None:
        """
        Adds all the given requirements to this execution parameter. Warns the
        user if incompatible requirements are detected.
        :param new_requirements:
        :return:
        """
        for req in new_requirements:
            if req in self._requires and \
                    self._requires[req] != new_requirements[req]:
                s = f"Conflicting requirements detected." \
                    f"\n\tParams: {repr(self)}" \
                    f"\n\tOld Req: {self._requires[req]}" \
                    f"\n\tNew Req: {new_requirements[req]}"
                logger.error(s)
                print(s, file=sys.stderr)
            self._requires[req] = new_requirements[req]

    def get_requirements(self) -> Dict[str, str]:
        """
        :return: All the packages their versions required by these execution
        parameters.
        """
        return copy(self._requires)

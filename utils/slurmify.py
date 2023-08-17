import os
import sys
from datetime import timedelta, datetime
from io import TextIOWrapper
from pathlib import Path

import logging
from threading import Thread
from typing import List, Tuple, TextIO

from utils.job_formatter import Runtime, JobBuilder, Job
from utils.path_manager import PathManager, quotes

global paths_manager

logger = logging.getLogger('root')


class SlurmJob(Job):
    """
    === Description ===
    A class that allows for the execution of a slurm script that wraps a command rather that
    directly executing a command in the current shell.

    === Private Attributes ===
    cmd: The command that this job will run.
    slurm_script_cmd: A command that, when executed, will run <cmd> through a slurm script.
    """
    _slurm_script_cmd: str

    def __init__(self, cmd: str, slurm_script_cmd: str) -> None:
        super().__init__(cmd)
        self._slurm_script_cmd = slurm_script_cmd

    def execute(self) -> None:
        logger.info('Executing slurm script' + self._slurm_script_cmd)
        os.system(self._slurm_script_cmd)
        self.complete = True

class FileToConsole:
    """
    === Description ===
    Redirects the output from a text file to a pipe as it is being written.

    === Private Attributes ===
    _in_file: The filepath to read from.
    _out_file: The file number to write to.
    _reader_thread: The thread that actively reads from the target file.
    """
    _in_file: Path
    _out_file: TextIO
    _reader_thread: Thread

    def __init__(self, in_file: Path, out_file: TextIO) -> None:
        """
        Initializes a file to console redirector.
        :param infile: The filepath to read from.
        :param out: The file to write to.
        """
        self._in_file = in_file
        self._out_file = out_file
        self._reader_thread = Thread(target=self.reader_funciton)

    def reader_funciton(self) -> None:
        """
        Redirect input from <self._in_file> to <self._out_file>
        """
        with open(self._in_file, 'r') as file:
            while True:
                line = file.readline()
                print(line, self._out_file)

    def start(self) -> None:
        """
        Start redirecting output.
        """
        self._reader_thread.start()



class Slurmifier(JobBuilder):
    """
    === Description ===
    Run system commands using the SLURM job scheduler.

    === Public Attributes ===

    === Private Attributes ===
    _email: User email address.
    _mailtype: What job events should prompt an email to the user.
    _priority: Whether to add jobs to the priority queue.

    """
    _email: str
    _mailtype: str
    _priority: bool
    _path_manager: PathManager
    _slurm_scripts_dir: Path

    _std_out_reader: FileToConsole
    _std_err_reader: FileToConsole

    def __init__(self, path_manager: PathManager, mailtype: List[str] = None, email_address: str = '',
                 priority: bool = False, redirect_to_console: bool = True) -> None:

        if mailtype is None:
            self.mailtype = ['NONE']
        else:
            self.mailtype = mailtype

        if email_address is '':
            email_address = 'kadochslurm@gmail.com'

        self._email = email_address
        self._priority = priority
        self._path_manager = path_manager
        self._slurm_scripts_dir = self._path_manager.make_dir(Path('slurm_scripts'))

        self.sbatch_out_path = self._path_manager.project_dir / 'sbatch.out'
        self.sbatch_err_path = self._path_manager.project_dir / 'sbatch.err'

        self._std_out_reader = FileToConsole(self.sbatch_out_path, sys.stdout)
        self._std_err_reader = FileToConsole(self.sbatch_out_path, sys.stderr)

        if redirect_to_console:
            self._std_out_reader.start()
            self._std_err_reader.start()

    def get_rt_params(self, runtime: Runtime) -> Tuple[str, str]:
        """Get slurm parameters derived from the runtime.
        :param Runtime runtime: A tuple containing the maximal days, minutes, and hours for which the program is
        expected to run.

        :return Tuple: the queue in which the job should be placed, the runtime string"""

        runtime_string = f"{runtime.days:01}-{runtime.hours:02}:{runtime.minutes:02}"

        if self._priority:
            queue = 'priority'
        else:
            rt = timedelta(days=runtime.days, hours=runtime.hours, minutes=runtime.minutes)

            if rt < timedelta(hours=12):
                queue = 'short'
            elif rt < timedelta(days=5):
                queue = 'medium'
            else:
                queue = 'long'

        return queue, runtime_string

    def slurmify(self, command: str, num_cores: int, runtime: Runtime, memory_per_core: int) -> SlurmJob:
        """Create a slurm job using the given parameters.

        :param str command: The command to be executed.
        :param int num_cores: The number of cores to use for the job.
        :param Runtime runtime: The maximum runtime of the job.
        :param int memory_per_core: The memory to reserve per CPU, in MiB.
        :param str std_out: The filepath to direct program output to. Directs to stdout by default.
        :param str std_err: The filepath to direct program output to. Directs to stderr by default.

        :return Job: Returns an executable job."""

        queue, runtime_str = self.get_rt_params(runtime)

        prog = 'command'.split(' ')[0].split('/')[-1]

        slurm_script = f'''#!/bin/bash
        #SBATCH --job-name=PythonPipelineRunning{prog}
        #SBATCH --mail-type={','.join(self.mailtype)}
        #SBATCH --ntasks={num_cores}
        #SBATCH --mem-per-cpu={memory_per_core}M
        #SBATCH --mail-user={self._email}
        #SBATCH -p {queue}
        #SBATCH -t {runtime_str}
        #SBATCH -o {self.sbatch_out_path}
        #SBATCH -e {self.sbatch_err_path}
        
        {command}
        '''

        # Make a unique script name.
        slurm_script_name = datetime.now().strftime("%m-%d-%Y|%H-%M-%S.%f.sh")
        slurm_command = 'sbatch ' + quotes(self._slurm_scripts_dir / slurm_script_name)

        return SlurmJob(command, slurm_command)

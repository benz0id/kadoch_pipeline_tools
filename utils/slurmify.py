import io
import os
import sys
from contextlib import redirect_stdout
from datetime import timedelta, datetime
from io import TextIOWrapper, StringIO, BufferedIOBase, BytesIO, BufferedRWPair
from pathlib import Path

import logging
from threading import Thread, Event
from time import sleep
from typing import List, Tuple, TextIO, Any

from patterns.observer import Observer
from utils.job_formatter import Runtime, JobBuilder, Job, ExecParams
from utils.path_manager import PathManager, quotes

global paths_manager

# Core charge rate in dollars per hour.
CORE_CHARGE_RATE = 0.0105

# Cost of RAM in dollars per hour.
RAM_CHARGE_RATE = 0.0013

# How long to wait before trying to read from one of the input files.
READER_SLEEP_INTERVAL = 0.1

# How long to wait before checking if all active jobs have finished.
SLURMIFIER_WAIT_INTERVAL = 0.01

# Print out the jobs still running in the array every this interval.
PRINT_WAITING_FOR_EVERY = 5

logger = logging.getLogger(__name__)


def get_O2_cost(params: ExecParams) -> float:
    """
    Gets the maximum cost of running the job given by <params> on the HMS O2 cluster.
    :param params: The parameters of the job.
    :return: Maximum cost of running the job on O2 in USD.
    """
    hours = params.max_runtime.days * 24 + params.max_runtime.hours + params.max_runtime.minutes // 60
    cpu_cost = hours * params.num_cores * CORE_CHARGE_RATE
    ram_cost = hours * params.ram_per_core * params.num_cores * RAM_CHARGE_RATE / 1024
    return cpu_cost + ram_cost


class O2Job(Job):
    """
    === Description ===
    A class that allows for the execution of a slurm script that wraps a command rather that
    directly executing a command in the current shell.

    === Private Attributes ===
    cmd: The command that this job will run.
    slurm_script_cmd: A command that, when executed, will run <cmd> through a slurm script.

    === Public Attributes ===
    cost: This estimated maximum cost of running this job on the O2 Cluster.
    """
    _slurm_script_cmd: str
    cost: float

    def __init__(self, cmd: str, slurm_script_cmd: str, cost: float,
                 observers: List[Observer] = None) -> None:
        super().__init__(cmd, observers)
        self.cost = cost
        self._slurm_script_cmd = slurm_script_cmd

    def execute(self) -> None:
        logger.info('Executing slurm script' + self._slurm_script_cmd)
        os.system(self._slurm_script_cmd)
        self.complete = True
        self.notify_observers()


class ArrayableO2(O2Job):
    """
    === Description ===
    A variant O2 job that executes the given job in another thread in order to
    avoid pausing.

    === Public Attributes ===
    thread: the thread started by this job.
    """

    thread: Thread

    def __init__(self, cmd: str, slurm_script_cmd: str, cost: float,
                 observers: List[Observer] = None) -> None:
        super().__init__(cmd, slurm_script_cmd, cost, observers)

    def execute(self) -> None:
        logger.info('Executing slurm script' + self._slurm_script_cmd)
        self.thread = Thread(target=os.system, args=[self._slurm_script_cmd])
        self.thread.start()
        self.complete = False
        self.notify_observers()


class FileToConsole:
    """
    === Description ===
    Redirects the output from a text file to a pipe as it is being written.

    This is a fairly sketchy implementation - but it works.

    === Private Attributes ===
    _in_file: The filepath to read from.
    _out_file: The file number to write to.
    _reader_thread: The thread that actively reads from the target file.
    _stop_event: Used to terminate dependent threads.
    """
    _in_file: Path
    _out_file: TextIO
    _reader_thread: Thread
    _stop_event: Event
    _start_lines: int

    def __init__(self, in_file: Path, out_file: TextIO) -> None:
        """
        Initializes a file to console redirector.
        :param infile: The filepath to read from.
        :param out: The file to write to.
        """
        self._in_file = in_file
        self._out_file = out_file
        self._stop_event = Event()
        self._reader_thread = Thread(target=self.reader_function)
        with open(in_file, 'r') as file:
            self._start_lines = len(file.readlines())

    def reader_function(self) -> None:
        """
        Redirect input from <self._in_file> to <self._out_file>
        """
        with open(self._in_file, 'r') as file:
            while True:
                line = file.readline().strip()

                if self._start_lines > 0:
                    self._start_lines -= 1
                elif line:
                    print(line.strip(), file=self._out_file)
                elif self._stop_event.is_set():
                    break
                else:
                    sleep(READER_SLEEP_INTERVAL)

    def start(self) -> None:
        """
        Start redirecting output.
        """
        self._reader_thread.start()

    def stop(self) -> None:
        """
        Stop redirecting output.
        """
        self._stop_event.set()


class StdOutInterceptor:
    """
    Intercept data written to stdout, store it then re-output to stdout.

    If a line output to stdout is a sbatch job request, store the job ID.
    """
    intercepted_jobs: List[int]
    snoopy_mode_active: bool
    pos: int

    def __init__(self):
        self._stringio = StringIO()
        self._original_stdout = sys.stdout
        self._stop_event = Event()
        self._interceptor_thread = Thread(target=self.intercept)
        self.snoopy_mode_active = False
        self.pos = 0

    def intercept(self) -> None:
        while True:
            self._stringio.seek(self.pos)
            line = self._stringio.readline()
            if line:
                self.pos += len(line)
                line = line.strip()
                print(line, 'aaaaa', file=self._original_stdout)
            if line and "Submitted batch job " in line \
                    and self.snoopy_mode_active:
                job_id = line.strip().split(' ')[-1]
                self.intercepted_jobs.append(int(job_id))
            elif self._stop_event.is_set():
                return
            else:
                sleep(READER_SLEEP_INTERVAL)

    def start(self) -> None:
        """
        Start redirecting output.
        """
        redirect_stdout(self._stringio)
        self._interceptor_thread.start()

    def stop(self) -> None:
        """
        Stop redirecting output.
        """
        self._stop_event.set()
        redirect_stdout(self._original_stdout)
        self._stop_event = Event()
        self._interceptor_thread = Thread(target=self.intercept)


class Slurmifier(JobBuilder, Observer):
    """
    === Description ===
    Run system commands using the SLURM job scheduler.

    === Private Attributes ===

    email: User email address.
    mailtype: What job events should prompt an email to the user.
    priority: Whether to add jobs to the priority queue.
    path_manager: Source for project path information.
    slurm_scripts_dir: Directory containing temporary slurm scripts created
        by this Slurmifier.

    active_threads: A list of threads stated by this Slurmifier, each
        executing some command. Used to enable waiting for the completion of
        job arrays.

    array_mode_active: If true, all slurm jobs are executed in seperate
        threads, and are made to pause. This allows for running and
        progress tracking of multiple slurm jobs.

    std_out_reader: Reads input to <sbatch_out_path> to stdout if active.
    std_err_reader: Reads input to <sbatch_err_path> to stderr if active.

    max_cost: The maximum aggregate cost of all jobs pushed to this class.
    max_current_expenditure: The maximum amount that could have possibly been
        spent by this slurmifier throughout it's execution. This is likely a
        wicked overestimate.


    === Public Attributes ===

    sbatch_out_path: Path to the sbatch output text file.
    sbatch_err_path: Path to the sbatch error text file.
    """
    _email: str
    _mailtype: str
    _priority: bool
    _path_manager: PathManager
    _slurm_scripts_dir: Path

    _active_threads: List[Thread]
    _array_mode_active: bool

    sbatch_out_path: Path
    sbatch_err_path: Path

    _std_out_reader: FileToConsole
    _std_err_reader: FileToConsole

    _max_cost: float
    _max_current_expenditure: float

    def __init__(self, path_manager: PathManager, mailtype: List[str] = None,
                 email_address: str = '', priority: bool = False,
                 redirect_to_console: bool = True, max_cost: float = 20
                 ) -> None:
        if mailtype is None:
            self.mailtype = ['NONE']
        else:
            self.mailtype = mailtype

        if email_address == '':
            email_address = 'kadochslurm@gmail.com'

        self._email = email_address
        self._priority = priority
        self._path_manager = path_manager
        self._slurm_scripts_dir = self._path_manager.rel_make(
            Path('.slurm_scripts'))

        self._active_threads = []
        self._array_mode_active = False

        self._max_cost = max_cost
        self._max_current_expenditure = 0

        self.config_writefiles()

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
            rt = timedelta(days=runtime.days, hours=runtime.hours,
                           minutes=runtime.minutes)

            if rt < timedelta(hours=12):
                queue = 'short'
            elif rt < timedelta(days=5):
                queue = 'medium'
            else:
                queue = 'long'

        return queue, runtime_string

    def config_writefiles(self) -> None:
        """Configure the slurm out and err files."""
        self.sbatch_out_path = self._path_manager.project_dir / 'sbatch.out'
        self.sbatch_err_path = self._path_manager.project_dir / 'sbatch.err'

        if not self.sbatch_out_path.exists() or not self.sbatch_err_path.exists():
            self.sbatch_out_path.touch()
            self.sbatch_err_path.touch()

            with open(self.sbatch_out_path, 'w') as out, open(
                    self.sbatch_out_path, 'w') as err:
                out.write('\n')
                err.write('\n')

    def slurmify(self, command: str, params: ExecParams) -> O2Job:
        """Create a slurm job using the given parameters.

        :param str command: The command to be executed.
        :param params: The parameters used to configure the environment in which the
            job is executed.

        :return Job: Returns an executable job."""

        queue, runtime_str = self.get_rt_params(params.max_runtime)

        prog = command.split(' ')[0].split('/')[-1]

        slurm_script = \
            [f'#!/bin/bash',
             f'#SBATCH --job-name={prog}',
             f'#SBATCH --mail-type={",".join(self.mailtype)}',
             f'#SBATCH --ntasks={params.num_cores}',
             f'#SBATCH --mem-per-cpu={params.ram_per_core}M',
             f'#SBATCH --mail-user={self._email}',
             f'#SBATCH -p {queue}',
             f'#SBATCH -t {runtime_str}',
             f'#SBATCH -o {self.sbatch_out_path}',
             f'#SBATCH -e {self.sbatch_err_path}'
             ]

        if params.wait or self._array_mode_active:
            slurm_script.append(f'#SBATCH -W')

        # Add lines to install the given modules.
        for requirement in params.get_requirements():
            slurm_script.append(
                f'module load {requirement}/{params.get_requirements()[requirement]}')

        slurm_script.append(command)

        slurm_script = '\n'.join(slurm_script)

        # Make a unique script name.
        slurm_script_name = datetime.now().strftime("%m-%d-%Y|%H-%M-%S.%f.sh")
        with open(self._slurm_scripts_dir / slurm_script_name, 'w') as file:
            file.write(slurm_script)
        slurm_command = 'sbatch ' + quotes(
            self._slurm_scripts_dir / slurm_script_name)

        logger.info("Slurm Script Successfully Created")
        logger.info(command + ' -> ' + slurm_command)
        logger.debug(slurm_script)

        if self._array_mode_active:
            return ArrayableO2(command, slurm_command, get_O2_cost(params),
                               observers=[self])
        else:
            return O2Job(command, slurm_command, get_O2_cost(params),
                         observers=[self])

    def notify(self, obj: Any) -> None:
        """
        === Notify Behaviours ===
        If is_instance(O2Job, obj):
            Job has been executed. Add cost to net cost.

        """
        if isinstance(obj, ArrayableO2):
            self._active_threads.append(obj.thread)
        if isinstance(obj, O2Job):
            logger.info(f'Updating current expenditure from '
                        f'{self._max_current_expenditure} to '
                        f'{self._max_current_expenditure + obj.cost}')
            self._max_current_expenditure += obj.cost

    def begin_array(self) -> None:
        """
        All subsequent slurm jobs are started in separate threads and are made
        to be pausing. Call <self.wait_for_all_jobs> in order to wait for all
        of these jobs to complete.
        """
        self._array_mode_active = True

    def wait_for_all_jobs(self) -> None:
        """
        Wait for all slurm jobs started by the slurmifier since the last
        <self.begin_array> call to complete.
        """

        if not any([thread.is_alive() for thread in self._active_threads]):
            self._array_mode_active = False
            return
        t = PRINT_WAITING_FOR_EVERY

        alive = []
        for thread in self._active_threads:
            if thread.is_alive():
                alive.append(thread)
        while alive:
            sleep(SLURMIFIER_WAIT_INTERVAL)
            t -= SLURMIFIER_WAIT_INTERVAL

            completed = []

            for thread in alive:
                if not thread.is_alive():
                    alive.remove(thread)
                    completed.append(thread)

            if completed:
                for thread in completed:
                    print('Complete: ',  str(thread))

        print('All threads complete.')

        self._array_mode_active = False

    def drop_array(self) -> None:
        """
        Disable array mode, and delete all active jobs.
        :return:
        """
        self._array_mode_active = False
        self._active_threads = []

    def prepare_job(self, cmd: str, exec_params: ExecParams) -> O2Job:
        """
        Creates a job that will execute the given <cmd> in an environment
        that provides the resources as specified in <exec_params>.
        :param cmd: The command to be executed.
        :param exec_params: Parameters for the runtime environment of the
        command.
        :return: An executable job.
        """
        self.do_expenditure_check(cmd, exec_params)

        job = self.slurmify(cmd, exec_params)
        logger.info('Slurm job successfully created.')
        return job

    def do_expenditure_check(self, cmd: str, params: ExecParams) -> None:
        """
        Checks that executing a job with <params> would not push this slurmifier beyond it's cost parameters.
        If so, stop the program.
        :param cmd: The command to be executed.
        :param params: The job parameters to be checked.
        """
        est_max_cost = get_O2_cost(params)
        logger.info(f"\n\t\t              Max Job cost: {est_max_cost:6.2f}\n"
                    f"\t\t       Current Expenditure: {self._max_current_expenditure:6.2f}\n"
                    f"\t\tRemaining after Completion: {self._max_cost - self._max_current_expenditure - est_max_cost:6.2f}")

        if self._max_current_expenditure + est_max_cost > self._max_cost:
            overage = -1 * self._max_cost - self._max_current_expenditure - est_max_cost
            logger.error(
                "Executing job would exceed allocated spending amount for this scheduler.")
            print(
                'Received request for job that would exceed allocated O2 spending.' +
                f"\n\n\t Command: {cmd}\n\n"
                f"\t\t                Max Cost of Job: {est_max_cost:6.2f}\n"
                f"\t\t                 Allotted Funds: {self._max_cost:6.2f}\n"
                f"\t\t            Current Expenditure: {self._max_current_expenditure:6.2f}\n"
                f"\t\t Overage that Would be Incurred: {overage:6.2f}"
                , file=sys.stderr)
            exit(1)

    def expenditure_report(self) -> None:
        """
        Print a quick report of the cost of the jobs run by this manager.
        """
        print(f"Expenditure - max ${self._max_current_expenditure:.2f} / "
              f"$ {self._max_cost:.2f}")

    def stop_readers(self) -> None:
        """
        Stops the threads reading from the slurm output files.
        """
        self._std_out_reader.stop()
        self._std_err_reader.stop()

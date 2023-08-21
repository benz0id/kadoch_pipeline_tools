import json
import logging
import os
from pathlib import Path
from typing import Union, List

from paramiko.client import SSHClient
from paramiko.sftp_client import SFTPClient

from utils.job_formatter import PythonJob
from utils.path_manager import quotes

logger = logging.getLogger(__name__)

# https://stackoverflow.com/questions/1392413/calculating-a-directorys-size-using-python
def get_size(start_path: str) -> int:
    """
    Get the size of a local directory.
    :param start_path: Directory containing all subdirectories and files to be sized.
    :return: The size of the directory at <start_path>.
    """
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            # skip if it is symbolic link
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)

    return total_size


class SSHFileTransferManager:
    """
    === Description ===

    Fetches and sends files to and from external machines using the SFTP through an ssh session that can be kept open,
    requiring minimal login requests.

    Extends and combines functionality from a paramiko.SSHClient and paramiko.SFTPClient

    === Private Attributes ===
    client: The active ssh client.
    hostname: The hostname of the client to which this manager is connected.
    sftp: The SFTP client used to tranfer files to and from the host.
    """
    _hostname: str
    _client: SSHClient
    _sftp: SFTPClient

    def __init__(self, config_json: Path) -> None:
        """
        Configured using a json file that may have the following format:

        {
            "host": "<hostname>",
            "username": "<username>",
            "password": "<password>"
        }

        Uses local .ssh key for login too.

        :param config_json: Path to the json file to be used.
        """
        self._client = SSHClient()
        self._client.load_system_host_keys()
        self._client.set_log_channel(__name__)
        self.connect_via_json(config_json)
        self._sftp = self._client.open_sftp()

    def connect_via_json(self, json_path: Path) -> None:
        """
        Attempts to connect to the server using credentials stored in the given json file.
        the json file must contain parameters "username", "password", and "host.

        :param json_path: Path to json file containing credentials.
        """
        file = open(json_path, 'r')
        data = json.load(file)
        self._hostname = data['host']
        self._client.connect(hostname=data['host'],
                             username=data['username'],
                             password=data['password'])

    def get_remote_size(self, filepath: Path) -> int:
        """
        Gets the size of the file or directory at <filepath> on the remote server.
        :param filepath: The path to the file or directory.
        :return:
        """

        if self.remote_is_file(filepath):
            stdin, stdout, stderr = self._client.exec_command("wc -c " + str(filepath))
            out = stdout.readline().strip()
            size = out.split(' ')[0]
            size = size.split('\t')[0]
            return int(size)
        elif self.remote_is_dir(filepath):
            stdin, stdout, stderr = self._client.exec_command("du -s -B1 " + str(filepath))
            out = stdout.readline().strip()
            size = out.split(' ')[0]
            size = size.split('\t')[0]
            return int(size)

    def get_local_size(self, local_filepath: Path):
        if local_filepath.is_dir():
            return get_size(start_path=str(local_filepath))
        else:
            return os.path.getsize(local_filepath)

    def remote_is_dir(self, filepath):
        cmd = f"[ -d {quotes(filepath)} ] && echo Exists"
        stdin, stdout, stderr = self._client.exec_command(cmd)
        out = stdout.readline().strip()
        is_dir = 'Exists' in out
        return is_dir

    def remote_is_file(self, filepath):
        cmd = f"[ -f {quotes(filepath)} ] && echo Exists"
        stdin, stdout, stderr = self._client.exec_command(cmd)
        out = stdout.readline().strip()
        is_dir = 'Exists' in out
        return is_dir

    def get_files_ignore_existing(self, remote_files: Union[Path, List[Path]],
                                  local_dir: Path) -> None:
        """
        Fetches all the given directories/files in <remote_files> and places the into <local_dir>. Files
        that already exist in <local_dir> are ignored.
        :param remote_files: The files to be copied locally.
        :param local_dir: The directory that will contain the copied files.
        :return:
        """

        if not isinstance(remote_files, list):
            remote_files = [remote_files]

        to_transfer = []

        for remote_filepath in remote_files:
            filename = str(remote_filepath).split('/')[-1]

            logger.info(f'Evaluating whether to fetch {filename} from {self._hostname}.')

            exists = (local_dir / filename).exists()

            remote_size = self.get_remote_size(remote_filepath)
            if exists:
                local_size = self.get_local_size(local_dir / filename)
            else:
                local_size = 0

            logger.info(f'\n\t\t Exists locally: {exists} \n\t\t Local size: {local_size} '
                        f'\n\t\t Remote size: {remote_size}')

            if not exists or local_size < remote_size:
                logger.info("Transfer will proceed.")
                to_transfer.append(remote_filepath)
            else:
                logger.info("Transfer will not proceed.")

        for remote_filepath in to_transfer:
            remote_filename = str(remote_filepath).split('/')[-1]
            if self.remote_is_file(remote_filepath):
                self._sftp.get(str(remote_filepath), str(local_dir / remote_filename))
            else:
                self.get_dir(remote_filepath, local_dir)

    def get_files_ignore_existing_job(self, remote_files: Union[Path, List[Path]],
                                      local_dir: Path) -> PythonJob:
        return PythonJob(cmd=f'get_files {str(remote_files)} -> {str(local_dir)}',
                         observers=[],
                         to_execute=self.get_files_ignore_existing,
                         remote_files=remote_files,
                         local_dir=local_dir)

    def put_files_ignore_existing_job(self, local_files: Union[Path, List[Path]],
                                      remote_dir: Path) -> PythonJob:
        return PythonJob(cmd=f'get_files {str(local_files)} -> {str(remote_dir)}',
                         observers=[],
                         to_execute=self.put_files_ignore_existing,
                         local_files=local_files,
                         remote_dir=remote_dir)

    def put_files_ignore_existing(self, local_files: Union[Path, List[Path]],
                                  remote_dir: Path) -> None:
        """
        Fetches all the given directories/files in <remote_files> and places the into <local_dir>. Files
        that already exist in <local_dir> are ignored.
        :param local_files: The files to be copied remotely.
        :param remote_dir: The remote directories into which those files should be copied.
        :return:
        """

        if not isinstance(local_files, list):
            local_files = [local_files]

        to_transfer = []

        for local_filepath in local_files:
            filename = str(local_filepath).split('/')[-1]

            logger.info(f'Evaluating whether to send {filename} to {self._hostname}.')

            hypothetical_path = remote_dir / filename
            exists = self.remote_is_dir(hypothetical_path) or self.remote_is_file(hypothetical_path)

            local_size = self.get_local_size(local_filepath)

            if exists:
                remote_size = self.get_remote_size(hypothetical_path)
            else:
                remote_size = 0

            logger.info(f'\n\t\t Exists locally: {exists} \n\t\t Local size: {local_size} '
                        f'\n\t\t Remote size {remote_size}')

            if not exists or local_size > remote_size:
                logger.info("Transfer will proceed.")
                to_transfer.append(local_filepath)
            else:
                logger.info("Transfer will not proceed.")

        for local_filepath in to_transfer:
            local_filename = str(local_filepath).split('/')[-1]
            logger.info(f'Transferring {local_filepath} -> {remote_dir / local_filename}')
            if not local_filepath.is_dir():
                self._sftp.get(str(local_filepath), str(remote_dir / local_filename))
            else:
                self.put_dir(local_filepath, remote_dir)

    def put_dir(self, local_dir: Path, remote_dir: Path) -> None:
        """
        Copy <local_dir> as a subdirectory <remote_dir>.
        :param local_dir: The local directory to be moved.
        :param remote_dir: The remote folder to contain the copy of <local_dir>.
        :return:
        """
        local_dir_name = str(local_dir).split('/')[-1]
        self.mkdir(str(remote_dir / local_dir_name), ignore_existing=True)
        for item in os.listdir(local_dir):
            if os.path.isfile(local_dir / item):
                self._sftp.put(str(local_dir / item), str(remote_dir / local_dir_name / item))
            else:
                self.mkdir(str(remote_dir / item), ignore_existing=True)
                self.put_dir(local_dir / item, remote_dir / local_dir_name / item)

    def get_dir(self, remote_dir: Path, local_dir: Path) -> None:
        """
        Copy <remote_dir> as a subdirectory of <local_dir>
        :param local_dir: The local directory to be moved.
        :param remote_dir: The remote folder to contain the copy of <local_dir>.
        :return:
        """
        remote_dir_name = str(remote_dir).split('/')[-1]
        if not (local_dir / remote_dir_name).exists() or (local_dir / remote_dir_name).is_file():
            os.makedirs(str(local_dir / remote_dir_name), exist_ok=True)
        for item in self._sftp.listdir_attr(str(remote_dir)):
            filename = item.filename
            if self.remote_is_file(remote_dir / filename):
                self._sftp.get(str(remote_dir / filename), str(local_dir / remote_dir_name / filename))
            else:
                os.mkdir(local_dir / filename)
                self.get_dir(remote_dir / filename, local_dir / remote_dir_name / filename)

    def mkdir(self, path, mode=511, ignore_existing=False):
        """
        Make a directory at path, do nothing if that directory already exists iff <ignore_existing>.
        :param path: Path to the directory to be created.
        :param mode: Mode with which to make the directory :)
        :param ignore_existing: Whether to ignore existing repositories.
        :return:
        """
        try:
            self._sftp.mkdir(path, mode)
        except IOError:
            if ignore_existing:
                pass
            else:
                raise

    def close(self) -> None:
        self._client.close()
        self._sftp.close()


import json
from pathlib import Path
from typing import Union

from paramiko.client import SSHClient
from paramiko.ssh_exception import SSHException

from utils.path_manager import PathManager


class SSHFileTransferManager:
    """
    === Description ===

    Fetches files from other servers using the SFTP through an ssh session that can be kept open,
    requiring minimal login requests.

    Adapter class for paramiko.SSHClient

    === Private Attributes ===
    path_manager: Manages paths required for transfers.
    client: The active ssh client.
    """

    def __init__(self, pathmanager: PathManager, config_json: Path) -> None:
        self._client = SSHClient()
        self._client.load_system_host_keys()
        self._client.set_log_channel(__name__)
        self.connect_via_json(config_json)
        self.sftp = self._client.open_sftp()

    def connect_interactive(self) -> None:
        """
        Fetch connection username and password from the user and use them to connect to the server.
        :return:
        """
        retry_str = ''
        success = False

        while retry_str == '':
            username = input("Enter your username:")
            password = input("Enter your password:")

            success = self.safe_login(username, password)
            if not success:
                retry_str = input('Connection failed due to authentication error. Press enter to retry.').strip()

        if not success:
            print('"Continuing without validated login.')
        else:
            print("Connection successfully established.")

    def connect_via_json(self, json_path: Path) -> None:
        """
        Attempts to connect to the server using credentials stored in the given json file.
        the json file must contain parameters "username", "password", and "host.

        :param json_path: Path to json file containing credentials.
        """
        file = open(json_path, 'r')
        data = json.load(file)
        self._client.connect(hostname=data['host'],
                             username=data['username'],
                             password=data['password'])

    def safe_login(self, username: Union[str, None] = None,
                   password: Union[str, None] = None) -> bool:
        """
        Attempts to connect to <self._remote_host> using the given authentication information.
        :param username: Username to connect with.
        :param password: Password to connect with.
        :return: Whether the login failed due to an authetication error.
        """
        try:
            self._client.connect(self._remote_host, username=username, password=password)
        except SSHException:
            return False
        return True






import abc
from typing import Dict


class ProgramWrapper(abc.ABC):
    """
    A class intended to provide a python interface for an external program.
    """

    @abc.abstractmethod
    def _get_dependencies(self) -> Dict[str, str]:
        """Returns all the packages/modules/softwares that must be present
        in the environment in order for the program wrapped by this class to
        execute."""
        pass

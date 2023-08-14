from abc import ABC, abstractmethod


class RerunAvoidant(ABC):
    @abstractmethod
    def outputs_exist(self):
        pass

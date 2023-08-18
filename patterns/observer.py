from abc import ABC, abstractmethod
from typing import List, Any


class Observer(ABC):

    @abstractmethod
    def notify(self, obj: Any):
        pass


class Observable(ABC):

    _observers: List[Observer]

    def __init__(self, observers: List[Observer]) -> None:
        if observers is None:
            observers = []

        self._observers = observers[:]

    def add_observer(self, observer: Observer) -> None:
        self._observers.append(observer)

    def notify_observers(self) -> None:
        for observer in self._observers:
            observer.notify()


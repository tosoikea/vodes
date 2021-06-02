from abc import ABC, abstractmethod

class Analysis(ABC):
    @abstractmethod
    def __init__(self, problem):
        self._problem = problem

    @abstractmethod
    def absolute(self):
        pass
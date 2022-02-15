# from sage.all import *

class ObjectWithPartitions:
    def __init__(self,k):
        self._k = k
        self._s = Partitions(self.k).cardinality()
    @property
    def s(self):
        return self._s
    @s.setter
    def s(self, value):
        raise AttributeError('The attribute s cannot be re-assigned')

    def number_of_partitions(self):
        return self.s

    @property
    def k(self):
        return self._k
    @k.setter
    def k(self, value):
        raise AttributeError('The attribute k cannot be re-assigned')
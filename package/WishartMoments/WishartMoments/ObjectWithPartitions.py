from sage.all import *

class ObjectWithPartitions:
    def __init__(self,k):
        self._k = k
        self._n = Partitions(self.k).cardinality()
    @property
    def n(self):
        return self._n
    @n.setter
    def n(self, value):
        raise AttributeError('The attribute n cannot be re-assigned')

    def number_of_partitions(self):
        return self.n

    @property
    def k(self):
        return self._k
    @k.setter
    def k(self, value):
        raise AttributeError('The attribute k cannot be re-assigned')
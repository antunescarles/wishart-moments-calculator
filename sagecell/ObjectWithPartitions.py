# Copyright 2022 Carlos S. Antunes Percincula
# This file is part of WishartMoments.
#
# WishartMoments is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# WishartMoments is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with WishartMoments. If not, see <https://www.gnu.org/licenses/>. 

#from sage.all import *

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
# Copyright 2022 Carlos S. Antunes Percincula
# This file is part of WishartMoments.
#
# WishartMoments is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# WishartMoments is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with WishartMoments. If not, see <https://www.gnu.org/licenses/>. 

#from sage.all import *

class Jacks:
    def __init__(self,k):
        self.z = SymmetricFunctions(QQ).zonal()
        self.p = SymmetricFunctions(QQ).power()
        self.m = SymmetricFunctions(QQ).monomial()

        self.k = k
        self.s = Partitions(k).cardinality()

        self.P = Partitions(self.k).list()
#         self.P.reverse()

    def jack_polynomial(self,t):
        jackM = self.m(self.z[self.P[t]])


        coefm1 = jackM.coefficients()[0] # coefficient of m[1,..,1]

        jackM_Jnorm= (jackM*factorial(self.k))/coefm1 # J normalization of the Jack in the monomial basis

        jcoefs = {}
        jcoefs['m'] = jackM_Jnorm.coefficients()
        jcoefs['p'] = [self.p(jackM_Jnorm).coefficient(self.P[self.s-i-1]) for i in range(0,self.s) ]
        # Note: We should fix the method that builds the coef matrix to use the same order for t and for i

        return jcoefs
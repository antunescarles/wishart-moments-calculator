from sage.all import *

class Jacks2:
    def __init__(self,k):
        self.z = SymmetricFunctions(QQ).zonal()
        self.p = SymmetricFunctions(QQ).power()
        self.m = SymmetricFunctions(QQ).monomial()

        self.k = k
        self.n = Partitions(k).cardinality()

        self.P = Partitions(self.k).list()
#         self.P.reverse()

    def jack_polynomial(self,s):
        jackM = self.m(self.z[self.P[s]])


        coefm1 = jackM.coefficients()[0] # coefficient of m[1,..,1]

        jackM_Jnorm= (jackM*factorial(self.k))/coefm1 # J normalization of the Jack in the monomial basis

        jcoefs = {}
        jcoefs['m'] = jackM_Jnorm.coefficients()
        jcoefs['p'] = [self.p(jackM_Jnorm).coefficient(self.P[self.n-i-1]) for i in range(0,self.n) ]
        # Note: We should fix the method that builds the coef matrix to use the same order for s and for i

        return jcoefs
#from sage.all import *
#from .ObjectWithPartitions import *
#from .Jacks import *

import math # To use the method isnan() to check if variables are NaN or not.
import numpy as np
from bisect import bisect_left ###
from numpy.linalg import matrix_power

def decorator(self,*args):
    # We should ensure the expression that is passed is either a variable or a negative power of a variable. This won't work otherwise.
    matrix_var = args[0].variables()[0]

    matrix_power = args[1]

    # exponent of the matrix variable
    pair = args[0].coefficients(matrix_var)[0]

    if matrix_power == 1:
        matrix_part = latex(args[0])
    else:
        matrix_part = '%s^{%d}' % (latex(args[0]),  pair[1]*matrix_power)

    return '{(\\mathrm{tr} \\, %s)}' % (matrix_part)

def negative_exp_prettyfier(self,*args):
    return '{%s^{-%d}}' % (latex(args[0]), args[1])

def ev(self,*args):
    expr = args[0]
    matrix_var = args[0].variables()[0] # We should ensure the expr has only one variable
    pair = expr.coefficients(matrix_var)[0]

    # As the input can be something like (2*A)^(-1) we have to retrieve the negative exponent the var A already has,
    # and multiply it by the second argument.
    matrix_var_exponent = pair[1]

    return (pair[0]**args[1])*tr(matrix_var,matrix_var_exponent*args[1])

def negative_exp_prettyfier(self,*args):
    return '{%s^{-%d}}' % (latex(args[0]), args[1])

function('tr', print_latex_func = decorator , nargs = 2 )
function('trace', print_latex_func = decorator , nargs = 2 , eval_func = ev )

function('inv', print_latex_func = negative_exp_prettyfier, nargs = 2)

class Expectations(ObjectWithPartitions):

    # Dictionaries for substitution
    w = var('w')
    W = var('W')
    n = var('n',latex_name="n")
    S = var('S',latex_name="\\Sigma")
    Sinv = var('Sinv', latex_name = "\\Sigma")
    Winv = var('Winv',latex_name = "W")

    def __init__(self,k):
        super().__init__(k)

        self.P = Partitions(k).list()
        self.P.reverse() #distinguish this from the one with the reverse order somehow
        # We could add jacks as a instance variable

        # Rings we will use
        self.R0 = PolynomialRing(SR,'S')
        self.R0._latex_names = ['\\Sigma']
        
        self.R3 = PolynomialRing(QQ,'n,r')

        (v_L,v_L_inv, L, rr) = self.vectors_L()
        self.v_L = v_L
        
        self.v_L_inv = v_L_inv
        self.L = L
        self.rr = rr

        # Matrices B, D and M

        ## For any moment, that of W or W^-1 we will need Bk and its inverse.
        self.Bk = self.compute_Bk()
        self.IBk = self.Bk.inverse()

        ## We will compute Dk and Dstar_k* only if needed
        ## and store them in a dictionary with only two keys: '+' and '-' that wich values will be D and Dstar_k respectively.
        self.Dk = {}

        ## We will compute M = IBk*Dk*Bk and M* = Bk*(Dstar_k)*IBk only if needed
        ## and store them in a dictionary with only two keys: '+' and '-' that wich values will be M and M* respectively.
        self.M = {}

        # Dictionary for substitution in the right side of the equation (the one with Sigma)
        self.Dinv = {w : (2*S)**(-1)}

        self.Catalogue = {}
        self.CatalogueInv = {}

    def compute_Bk(self):
        t = 0 # partitions go like mu(s-(s-1)) = m(1) < mu(s-(s-2)) = mu(2) < mu(s-1) < mu(s-0) = [s]
        # So t can range from 0 to s-1
        # The program will compute the Jack polynomial corresponding to partition mu[s-t] of the list of partitions

        jacks = Jacks(self.k)

        ## Below we get a dictionary with 2 keys: 'p' for coeffs in power-sum basis, and 'm' for the coeffs in monomial basis.
        coef = jacks.jack_polynomial(t)
        Bk = matrix(QQ,1,coef['p']) # We use the ones of the power-sum basis.

        t=0
        t+=1
        while t <= self.s - 1: # we use while instead of for bc when k=2 range(1,1) is empty and it never enters the loop
            coef = jacks.jack_polynomial(t)
            row =  matrix(QQ,1,coef['p'])
            Bk = Bk.stack(row)
            t+=1

        return Bk

    def compute_Dk(self,inverse = False):

        P = Partitions(self.k).list()
        (n,r) = self.R3.gens()

        if not inverse:
            if not('+' in self.Dk):
                Dk = matrix(self.R3,self.s,self.s,0)

                pm = [1]*self.s
                for i in range(0,self.s):
                    lm = len(P[i])
                    for j in range(1,lm+1):
                            for t in range(1,P[i][j-1]+1):
                                pm[i] *= (n/2) +t-1- (j-1)/2
                    Dk[i,i] = pm[i]

                self.Dk['+'] = Dk
        else:
            if not('-' in self.Dk):
                ## Compute the diagonal for the expectations of the inverse
                R3_frac = self.R3.fraction_field()
                Dk_star = matrix(R3_frac,self.s,self.s,0)

                qm = [1]*self.s
                for i in range(0,self.s):
                    lm = len(P[i])
                    for j in range(self.k-lm+1,self.k+1):
                            for t in range(1,P[i][self.k-j+1 -1]+1):
                                qm[i] *= (n-r)/2 + (self.k-j+1)/2 -t 
                    denom = qm[i]
                    Dk_star[i,i] = 1/denom
                
                self.Dk['-'] = Dk_star

    def compute_M(self,inverse = False):

        self.compute_Dk(inverse)
        if not(inverse):
            if not('+' in self.M):
                self.M['+'] = self.IBk*self.Dk['+']*self.Bk
        else:
            if not('-' in self.M):
                self.M['-'] = self.IBk*self.Dk['-']*self.Bk

    def prettify_negative_powers_of_matrix_var(self, expr, matrix_var):
        ## Artifact to print E[\Sigma ^{-1}] nicely (if we don't do this Sigma^{-1} is printed as 1/Sigma which isn't pretty for a matrix)
        # 1) Extract the coefficients of every negative power of Sinv
        # 2) Form a new expression multiplying the coef of the (-j)-th powe of Sinv for a new variable, something like Sj with latex_name \Sigma^{-j}

        pairs = expr.coefficients(matrix_var) # we get the list of pairs of the form (coefficient of Sinv^power, power)

        # To do: change the name of the variable added.
        # Use Sinvj for Sinv^(-j) instead of S, and it will be probably needed to change the trace_decorator_inv
        expr = sum( [ p[0].factor()*inv(matrix_var,abs(p[1])) for p in pairs] ) # factorize the denominator

        return expr

    def moment(self, t, inverse = False):

        self.compute_M(inverse) # Computes M['+'] or M['-'] only if it hasn't already been computed.
        
        portrait = self.partition_to_portrait(self.P[t])

        if inverse :
            if self.P[t] in self.CatalogueInv :
                m = self.CatalogueInv[self.P[t]]
            else :
                variable = (self.compute_L(portrait,False)/self.k).subs({w:self.W**(-1)}).substitute_function(tr,trace)
                variable = self.prettify_negative_powers_of_matrix_var(variable, W)
                
                expectation = ((self.M['-'].row(t)*self.v_L_inv)/ self.k).subs(self.Dinv).substitute_function(tr,trace)
                expectation = self.prettify_negative_powers_of_matrix_var(expectation,S)
                
                m = { 'var':variable , 'moment':expectation }
                self.CatalogueInv[self.P[t]] = m

            return m

        if self.P[t] in self.Catalogue :
            m = self.Catalogue[self.P[t]]
        else :
            variable = (self.compute_L(portrait,False)/self.k).subs({w:W}).substitute_function(tr,trace)

            L = self.s*[NaN]
            for j in range(0,self.s):
                portrait = self.partition_to_portrait(self.P[j])
                L[j] = self.compute_L(portrait,True)
            
            expectation = sum([ self.M['+'][t,j]*L[j] for j in range(0,self.s)])/self.k

            m = { 'var':variable , 'moment':expectation }
            
            self.Catalogue[self.P[t]] = m

            
        return m

    def expressions(self,inverse=False):
        if inverse:
            v = (self.v_L_inv/self.k).subs({w:self.W**(-1)})

            var_list = []
            for i in range(0,self.s):
                  print([i,self.prettify_negative_powers_of_matrix_var(v[i].substitute_function(tr,trace), W)])
        else:
            v = (self.v_L/self.k).subs({w:W})
            var_list = []
            for i in range(0,self.s):
                print([i,v[i].substitute_function(tr,trace)])
                
    def expression(self, t,inverse=False):
        expr = []
        if not(inverse):
            expr = (self.v_L[t]/self.k).subs({w:W}).substitute_function(tr,trace)
        else:
            expr = (self.v_L_inv[t]/self.k).subs({w:self.W**(-1)})
            expr = self.prettify_negative_powers_of_matrix_var(expr.substitute_function(tr,trace), W)
        return expr
    def partition_to_portrait(self, lam):
        #  lam (lam) is a partition
        lam = list(lam) # we have to ensure we work with a list and not a object of another data type.

        i = [0]*self.k # i will represent the portrait (i) associated to lambda
        set_lam = set(lam)
        for j in set_lam:
            #  we want to represent to store st such that st[0]*1 + st[1]*2 + st[2]*3 + ... + st[k-1]*k = k
            # notice that index starts from zero but is the same. That's the reason why we add 1 to i below:
            i[j-1] = list(lam).count(j)
        return i

    def trace_decorator(self, l, varname):
        # l will be j+1, the power of the argument
        # p will be i[j], the power of the trace
        a = "\\mathrm{tr}\\,"
        if (l == 1):
            a = a + varname
        else:
            a = a+varname+"^%d"%(l)

        return "("+a+")"

    def trace_decorator_inv(self, l, varname):
        # l will be j+1, the power of the argument
        # p will be i[j], the power of the trace
        a = "\\mathrm{tr}\\,"+varname+"^{-%d}"%(l)

        return "("+a+")"

    def compute_r(self, i, right=False):
        # i is a portrait (i1, i2,..., ik)

        if not(right): # in this case we compute r_i for the lefthand-side of the equation
            w = var('w')
            
            r_i = prod([trace( w , j+1)**(i[j]) for j in range(0,self.k) ])
            
        else: # in this case, we need the righthand-side of the equation, where we use 2*Sigma
            S = self.R0.gen()
            
            x = var('S',latex_name='\\Sigma')
            r_i = prod([trace( 2*x , j+1)**(i[j]) for j in range(0,self.k) ])
            
        return r_i

    def compute_L(self,i,right=False):
        if not(right):
            w = var('w')
        
            r_i = self.compute_r(i,False)

            L_i = sum([expand( r_i*(j+1)*i[j]*w**(j+1)/trace(w,j+1) ) for j in range(0,self.k) ])
            # ^  for some reason if we multiply r[i] outside sum([...]) simplification is not done right...
        else:
            S = self.R0.gen()    
            x = var('S',latex_name = '\\Sigma')
            
            r_i = self.compute_r(i,True)
            
            L_i = sum([expand( r_i*(j+1)*i[j]*(2*S)**(j+1)/trace(2*x,j+1) ) for j in range(0,self.k) ])
            # ^  for some reason if we multiply r[i] outside sum([...]) simplification is not done right...
        return L_i

    def compute_numerical_value_r(self,i,S):
        tr = [ np.trace(matrix_power(S,j+1)) for j in range(0,self.k)] 

        r_i = prod([ (tr[j])**(i[j]) for j in range(0,self.k) ])
        return (r_i , tr)

    def compute_numerical_value_L(self,i,S):

        (r_i,tr) = self.compute_numerical_value_r(i,S)

        L_i = sum([ r_i*(j+1)*i[j]*matrix_power(S,j+1)/tr[j] for j in range(0,self.k) ]) 
        # ^  for some reason if we multiply r[i] outside sum([...]) simplification is not done right...

        return L_i

    def vectors_L(self):

        rr = [] # this shouldn't be named 'r' bc it crashes with the name of the parameter r that represents de dimension of Sigma
        L = []
        for j in range(0,self.s):
            rr.append(self.partition_to_portrait(self.P[j]))
            L.append(self.compute_L(rr[j]))

        v_L = vector(SR,L)

        # The next one is the same than v_L, but for caution we use a fresh variable.
        # When we have to print them, we'll substitute the variable used in v_L_inv for sigma^-1 in the latex representation.
        v_L_inv = vector(SR,L)
        return (v_L , v_L_inv, L, rr)


    def numerical_L_vectors(self,Sigma):
        A = Sigma

        Lnum = [] # for the right side.
        L_inv_num = []
        for j in range(0,self.s):
            Lnum.append(self.compute_numerical_value_L(self.rr[j],2*A)) # For the right-side that is not symbolic.

            # For the right-side of the inverse that is not symbolic. We call the same function but with the inverse of A as parameter.
            L_inv_num.append(self.compute_numerical_value_L(self.rr[j],matrix_power(2*A ,-1)))

        return (Lnum,L_inv_num)

    def evaluated_wishart_expectations(self,Sigma,n_param,inverse):

        A = Sigma
        dim_Sigma = Sigma.shape[0]
        (n,r) = self.R3.gens()

        (Lnum,L_inv_num) = self.numerical_L_vectors(A)

        # For the righthandside we have carry the computations one by one because we cannot form a vector of matrices.
        Enum = [NaN]*self.s # Numerical (concrete) expectation

        # Concrete inverse
        # For the righthandside we have carry the computations one by one because we cannot form a vector of matrices.
        E_inv_num = [NaN]*self.s # Evaluated expectation
        
        if not(inverse):
            for i in range(0,self.s):
                Enum[i] = sum([self.M['+'][i,j].subs({n: n_param})*Lnum[j] for j in range(0,self.s)])
        else:
            if (n_param > 2*self.k + (dim_Sigma -1)): # Condition for the inverse to be calculated
                for i in range(0,self.s):
                    E_inv_num[i] = sum([self.M['-'][i,j].subs({n: n_param, r: dim_Sigma})*L_inv_num[j] for j in range(0,self.s)])

        return (Enum, E_inv_num)

    def substitute_with_W_inverse(self,Ik_indx):
        # Change Sinv^{-j} for Sj because it is best for pretty printing it
        expr_1 = self.v_L_inv[Ik_indx-1].subs({w : W**(-1)})/self.k # at this point this expression contains Winv as a variable.

        expr_1 = expr_1.substitute_function(tr,trace)

        expr_2 = self.E_inv[Ik_indx-1]/self.k

        random_variable_inv = latex(self.prettify_negative_powers_of_matrix_var(expr_1 , W))
        expected_value_inv = latex(self.prettify_negative_powers_of_matrix_var(expr_2 , S))

        return (random_variable_inv, expected_value_inv)

    def evaluate_moment(self, t, n_param, Sigma, inverse=False):
        self.compute_M(inverse) # Computes M['+'] or M['-'] only if it hasn't already been computed.

        A = Sigma
        dim_Sigma = Sigma.shape[0]

        # Avoid syntactic sugar!
        # R2.<f,p,r> = QQ['f,p,r']
#        self.R2 = PolynomialRing(QQ,'f,p,r')
#        (f,p,r) = self.R2.gens()

        (Enum, E_inv_num)= self.evaluated_wishart_expectations(Sigma, n_param,inverse)

        if inverse == False:
            variable = (self.v_L[t]/self.k).subs({w:W}).substitute_function(tr,trace)
            evaluated_expectation = Enum[t]
        else:
            variable = (self.v_L_inv[t]/self.k).subs({w:self.W**(-1)}).substitute_function(tr,trace)
            variable = self.prettify_negative_powers_of_matrix_var(variable, W)

            evaluated_expectation = E_inv_num[t]

        eval_m = { 'var': variable, 'moment': evaluated_expectation }
        
        return eval_m

    def pretty_print_eval_moment(self, t, n_param, Sigma, inverse = False):
        eval_m = self.evaluate_moment(t,n_param,Sigma,inverse)
        pretty_print(html(r'<p style="margin-top:2em; margin-bottom:2em; margin-left:4.5em">$ \mathbb{E}(%s) = %s $</p>' % (latex(eval_m['var']),latex(matrix(eval_m['moment']))) ))

    def number_of_expectations(self):
        return self.number_of_partitions()

    def pretty_print_moment(self,t,inverse=False):
        m = self.moment(t,inverse)
        pretty_print(html(r'<p style="margin-top:2em; margin-bottom:2em; margin-left:4.5em">$ \mathbb{E}(%s) = %s $</p>' % (latex(m['var']),latex(m['moment'])) ))
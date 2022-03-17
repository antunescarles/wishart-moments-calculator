from sage.all import *
from .ObjectWithPartitions import *
from .Jacks import *


import math # To use the method isnan() to check if variables are NaN or not.
import numpy as np
from bisect import bisect_left ###
from numpy.linalg import matrix_power

def decorator(self,*args):
    # We need to ensure that the expression that is passed is either a variable or a negative power of a variable. This won't work otherwise.
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
    matrix_var = args[0].variables()[0] # We must ensure the expr has only one variable
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
    N = var('n',latex_name="n")
    S = var('S',latex_name="\\Sigma")
    Sinv = var('Sinv', latex_name = "\\Sigma")
    Winv = var('Winv',latex_name = "W")

    def __init__(self,k):
        super().__init__(k)

        self.P = Partitions(k).list()
        self.P.reverse() #distinguish this from the one with the reverse order somehow
        # We could add jacks as a instance variable

        # Rings we will use
        self.R2 = PolynomialRing(QQ,'f,p,r')
        (f,p,r) = self.R2.gens()

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

#         (Mp, M_pnf_star) = self.symbolic_M_matrices()
#         self.Mp = Mp
#         self.M_pnf_star = M_pnf_star

        # Dictionaries for substitution in the right side of the equation (the one with Sigma)
        (f,p,r) = self.R2.gens()

        self.D = {p:n/2 , w:2*S , f:1/2 }
        self.Dinv = {p:n/2 , w : (2*S)**(-1)  , f: 1/2}

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
        (f,p,r) = self.R2.gens()

        if not inverse:
            if not('+' in self.Dk):
                Dk = matrix(self.R2,self.s,self.s,0)

                pm = [1]*self.s
                for i in range(0,self.s):
                    lm = len(P[i])
                    for j in range(1,lm+1):
                            for t in range(1,P[i][j-1]+1):
                                pm[i] *= p +t-1- (j-1)*f
                    Dk[i,i] = pm[i]

                self.Dk['+'] = Dk
        else:
            if not('-' in self.Dk):
                ## Compute the diagonal for the expectations of the inverse
                R2_frac = self.R2.fraction_field()
                Dk_star = matrix(R2_frac,self.s,self.s,0)

                n = var('n',latex_name="n")

                qm = [1]*self.s
                for i in range(0,self.s):
                    lm = len(P[i])
                    for j in range(self.k-lm+1,self.k+1):
                            for t in range(1,P[i][self.k-j+1 -1]+1):
                                qm[i] *= p + (self.k-j+1)*f -t # here I'd like to use another var, e.g, q instead of the same p,
                                                          # but as Ill inmediatelly substitute it's not worth the effort thinking a better solution.
                    # Evaluate the expr. in q = p-r*f
                    denom = (qm[i].subs({p : (p - r*f)})) # later we'll substitute for f = 1/2 (as f=1/2 is the value of f we're interested in)
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

    def symbolic_M_matrices(self):
        # Calculates M(p) and M^*(p,r,f) in term of the parameters p,r, and f, which will be subsituted later.

        outmost_verbose = False

        t = 0 # partitions go like mu(s-(s-1)) = m(1) < mu(s-(s-2)) = mu(2) < mu(s-1) < mu(s-0) = [s]
        # So t can range from 0 to s-1
        # The program will compute the Jack polynomial corresponding to partition mu[s-t] of the list mu of partitions

        jacks = Jacks(self.k)

        ## Below we get a dictionary with 2 keys: 'p' for coeffs in power-sum basis, and 'm' for the coeffs in monomial basis.
        coef = jacks.jack_polynomial(t)
        Bk = matrix(QQ,1,coef['p']) # We use the ones of the power-sum basis.

        t=0
        t+=1
        # for t in range(1,s-1):
        while t <= self.s - 1: # we use while instead of for bc when k=2 range(1,1) is empty and it never enters the loop
            coef = jacks.jack_polynomial(t)
            row =  matrix(QQ,1,coef['p'])
            Bk = Bk.stack(row)
            t+=1

        # Avoid sage's syntactic sugar
        #R2.<f,p,r> = QQ['f,p,r']
        R2 = PolynomialRing(QQ,'f,p,r')
        (f,p,r) = R2.gens()

        P = Partitions(self.k).list()

        Dk = matrix(R2,self.s,self.s,0)

        pm = [1]*self.s
        for i in range(0,self.s):
            lm = len(P[i])
            for j in range(1,lm+1):
                    for t in range(1,P[i][j-1]+1):
                        pm[i] *= p +t-1- (j-1)*f
            Dk[i,i] = pm[i]

            if outmost_verbose: print(P[i]," -->  ", pm[i])

        # Compute Mp
        IBk = Bk.inverse()

        Mp = IBk*Dk*Bk

        ## Compute the diagonal for the expectations of the inverse
        R2_frac = R2.fraction_field()
        Dk_star = matrix(R2_frac,self.s,self.s,0)

        if outmost_verbose:  print("Elementos de la diagonal de Dk factorizados\n")

        n = var('n',latex_name="n")

        qm = [1]*self.s
        for i in range(0,self.s):
            lm = len(P[i])
            for j in range(self.k-lm+1,self.k+1):
                    for t in range(1,P[i][self.k-j+1 -1]+1):
                        qm[i] *= p + (self.k-j+1)*f -t # here I'd like to use another var, e.g, q instead of the same p,
                                                  # but as Ill inmediatelly substitute it's not worth the effort thinking a better solution.
            # Evaluate the expr. in q = p-r*f
            denom = (qm[i].subs({p : (p - r*f)})) # later we'll substitute for f = 1/2 (as f=1/2 is the value of f we're interested in)
            Dk_star[i,i] = 1/denom

        # When it corresponds, compute M^*(p-rf) r = Partitions(k).cardinality() == s
        M_pnf_star = IBk*Dk_star*Bk

        ## Esto es para cheuqear un error nomas
        #         show("B_k = "+ latex(Bk))
        #         show("B_k^{-1} = "+ latex(IBk))

        #         DD = Dk_star.subs({p:n/2})
        #         pretty_print(html(r'<center>$D_k^*(\frac{s-r}{2}) = \begin{pmatrix}%s & %s \\ %s & %s \end{pmatrix}$</center>' % (latex(DD[0,0].factor()) , latex(DD[0,1]), latex(DD[1,0]) , latex(DD[1,1].factor()) ) ))

        return (Mp,M_pnf_star)

    def prettify_negative_powers_of_matrix_var(self, expr, matrix_var):
        ## Artifact to print E[\Sigma ^{-1}] nicely (if we don't do this Sigma^{-1} is printed as 1/Sigma which isn't pretty for a matrix)
        # 1) Extract the coefficients of every negative power of Sinv
        # 2) Form a new expression multiplying the coef of the (-j)-th powe of Sinv for a new variable, something like Sj with latex_name \Sigma^{-j}

        pairs = expr.coefficients(matrix_var) # we get the list of pairs of the form (coefficient of Sinv^power, power)

        # To do: change the name of the variable added.
        # Use Sinvj for Sinv^(-j) instead of S, and it will be probably needed to change the trace_decorator_inv
        expr2 = sum( [ p[0].factor()*inv(matrix_var,abs(p[1])) for p in pairs] ) # factorize the denominator

        return expr2

    def moment(self, t, inverse = False):

        self.compute_M(inverse) # Computes M['+'] or M['-'] only if it hasn't already been computed.

        if inverse :
            if self.P[t] in self.CatalogueInv :
                m = self.CatalogueInv[self.P[t]]
            else :
                variable = (self.v_L_inv[t]/self.k).subs({w:self.W**(-1)}).substitute_function(tr,trace)
                variable = self.prettify_negative_powers_of_matrix_var(variable, W)

                expectation = ((self.M['-'].row(t)*self.v_L_inv)/ self.k).subs(self.Dinv).substitute_function(tr,trace)
                expectation = self.prettify_negative_powers_of_matrix_var(expectation,S)
                self.CatalogueInv[self.P[t]] = (variable,expectation)

                m = (variable, expectation)
            return m

        if self.P[t] in self.Catalogue :
            m = self.Catalogue[self.P[t]]
        else :
            variable = (self.v_L[t]/self.k).subs({w:W}).substitute_function(tr,trace)
            expectation = ((self.M['+'].row(t)*self.v_L)/ self.k).subs(self.D).substitute_function(tr,trace)
            self.Catalogue[self.P[t]] = (variable,expectation)

            # m = (variable, expectation)
            m = { 'var':variable , 'moment':expectation }
        return m

    def expressions(self,inverse=False):
        if inverse:
            v = (self.v_L_inv/self.k).subs({w:self.W**(-1)})

            var_list = []
            for i in range(0,self.s):
#                 var_list.append([i,self.prettify_negative_powers_of_matrix_var(v[i].substitute_function(tr,trace), W)])
                  print([i,self.prettify_negative_powers_of_matrix_var(v[i].substitute_function(tr,trace), W)])
        else:
            v = (self.v_L/self.k).subs({w:W})
            var_list = []
            for i in range(0,self.s):
#                 var_list.append([i,v[i].substitute_function(tr,trace)])
                print([i,v[i].substitute_function(tr,trace)])
#         return var_list
    def expression(self, t,inverse=False):
        expr = []
        if not(inverse):
            expr = (self.v_L[t]/self.k).subs({w:W}).substitute_function(tr,trace)
        else:
            expr = (self.v_L_inv[t]/self.k).subs({w:self.W**(-1)})
            expr = self.prettify_negative_powers_of_matrix_var(expr.substitute_function(tr,trace), W)
        return expr
    def partition_to_portrait(self,t):
        #  t is a type, or equivalently, a partition
        t = list(t) # we have to ensure we work with a list and not a object of another data type.

        i = [0]*self.k
        set_t = set(t)
        for j in set_t:
            #  we want to represent to store st such that st[0]*1 + st[1]*2 + st[2]*3 + ... + st[k-1]*k = k
            # notice that index starts from zero but is the same. That's the reason why we add 1 to i below:
            i[j-1] = list(t).count(j)
        return i

    def trace_decorator(self, l, varname):
        # l sera j+1, la potencia del argumento
        # p sera i[j], la potencia de la traza
        a = "\\mathrm{tr}\\,"
        if (l == 1):
            a = a + varname
        else:
            a = a+varname+"^%d"%(l)

        return "("+a+")"

    def trace_decorator_inv(self, l, varname):
        # l sera j+1, la potencia del argumento
        # p sera i[j], la potencia de la traza
        a = "\\mathrm{tr}\\,"+varname+"^{-%d}"%(l)

        return "("+a+")"

    def compute_r(self, i):
        # i is a portrait

        w = var('w')

        # When we have b1 we want tr\sigma instead of tr(\sigma^1)
#         r_i = prod([var('b%d'%(j+1),latex_name = self.trace_decorator(j+1,"\\sigma") )^(i[j]) for j in range(0,self.k) ])

#         A = var('A')
        r_i = prod([trace( w , j+1)**(i[j]) for j in range(0,self.k) ])
        return r_i

    def compute_L(self,i):
        w = var('w')

        r_i = self.compute_r(i)

#         L_i = sum([expand( r_i*(j+1)*i[j]*w^(j+1)/var('b%d'%(j+1)) ) for j in range(0,self.k) ])

        L_i = sum([expand( r_i*(j+1)*i[j]*w**(j+1)/trace(w,j+1) ) for j in range(0,self.k) ])
        # ^  por alguna razon si multiplicamos r[i] afuera de sum([...]) no simplifica la bien la expresión...

        return L_i

    def compute_numerical_value_r(self,i,S):
        # tr = [ (S**(j+1)).trace() for j in range(0,self.k)]
        tr = [ np.trace(matrix_power(S,j+1)) for j in range(0,self.k)] # Numpy version


        r_i = prod([ (tr[j])**(i[j]) for j in range(0,self.k) ])
        return (r_i , tr)

    def compute_numerical_value_L(self,i,S):

        (r_i,tr) = self.compute_numerical_value_r(i,S)

        # L_i = sum([ r_i*(j+1)*i[j]*S**(j+1)/tr[j] for j in range(0,self.k) ])
        L_i = sum([ r_i*(j+1)*i[j]*matrix_power(S,j+1)/tr[j] for j in range(0,self.k) ]) # Numpy version

        # ^  por alguna razon si multiplicamos r[i] afuera de sum([...]) no simplifica la bien la expresión...

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

            # assert A.is_invertible() , "Error: A is not invertible." # think what happens if A is not over QQ, but over RR.

            # For the right-side of the inverse that is not symbolic. We call the same function but with the inverse of A as parameter.
            L_inv_num.append(self.compute_numerical_value_L(self.rr[j],matrix_power(2*A ,-1)))

        return (Lnum,L_inv_num)

    def wishart_expectations_numval(self,Sigma,n_param,inverse):

        A = Sigma
        dim_Sigma = Sigma.shape[0]

        outmost_verbose = False

        (Lnum,L_inv_num) = self.numerical_L_vectors(A)

        #R2.<f,p,r> = QQ['f,p,r']
        R2 = PolynomialRing(QQ,'f,p,r')
        (f,p,r) = R2.gens()

        #Para el lado derecho hay que hacer las cuentas mas a mano porque no podemos formar un vector de matrices...
        Enum = [NaN]*self.s # Numerical (concrete) expectation

        # Concrete inverse
        #Para el lado derecho hay que hacer las cuentas mas a mano porque no podemos formar un vector de matrices...
        E_inv_num = [NaN]*self.s # Numerical (concrete) expectation

        if not(inverse):
            for i in range(0,self.s):
                Enum[i] = sum([self.M['+'][i,j].subs({p: n_param/2, f:1/2 })*Lnum[j] for j in range(0,self.s)])
        else:
            if (n_param > 2*self.k + (dim_Sigma -1)): # Condition for the inverse to be calculated
                for i in range(0,self.s):
                    E_inv_num[i] = sum([self.M['-'][i,j].subs({p: n_param/2, r: dim_Sigma})*L_inv_num[j] for j in range(0,self.s)])

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
        self.R2 = PolynomialRing(QQ,'f,p,r')
        (f,p,r) = self.R2.gens()

        (Enum, E_inv_num)= self.wishart_expectations_numval(Sigma, n_param,inverse)

        if inverse == False:
            variable = (self.v_L[t]/self.k).subs({w:W}).substitute_function(tr,trace)
            # evaluated_expectation = Enum[t].subs({f:1/2})
            evaluated_expectation = Enum[t]
        else:
            variable = (self.v_L_inv[t]/self.k).subs({w:self.W**(-1)}).substitute_function(tr,trace)
            variable = self.prettify_negative_powers_of_matrix_var(variable, W)

            evaluated_expectation = E_inv_num[t].subs({f:1/2})

        # eval_m = (variable, evaluated_expectation)
        eval_m = { 'var': variable, 'moment': evaluated_expectation }

#         if self.P[t] in self.Catalogue :
#             m = self.Catalogue[self.P[t]]
#         else :
#             variable = (self.v_L[t]/self.k).subs({w:W}).substitute_function(tr,trace)
#             evaluated_expectation = Enum[t+1]
# #             expectation_value = ((self.M['+'].row(t)*self.v_L)/ self.k).subs(self.D).substitute_function(tr,trace)

#             eval_m = (variable, evaluated_expectation)
        return eval_m

    def pretty_print_eval_moment(self, t, n_param, Sigma, inverse = False):
        eval_m = self.evaluate_moment(t,n_param,Sigma,inverse)
        pretty_print(html(r'<p style="margin-top:2em; margin-bottom:2em; margin-left:4.5em">$ \mathbb{E}(%s) = %s $</p>' % (latex(eval_m[0]),latex(eval_m[1])) ))

#         (lsideD, new_E_inv_expr_lside, new_E_inv_expr) = self.expectations_expressions(Ik_indx)



#         pretty_print(html(r'<div>$(i) = %s $</div>' % LatexExpr(self.P[Ik_indx-1])) )
#         pretty_print(html(r'<p style= "margin-top:2em; margin-bottom:2em; margin-left:4.5em">$$2\Sigma = %s $$</p>' %latex(2*A) ))
#         pretty_print(html( r'<p style="margin-top:2em; margin-bottom:2em; margin-left:4.5em"> $$\mathbb{E}(%s) \; = \; %s$$</p>' % (latex(self.v_L[Ik_indx-1].subs(lsideD)/self.k) , latex(Enum[Ik_indx-1].subs({p:n/2})/self.k)) ))
#         pretty_print(html( r'$\text{The moments of } W^{-1} \text{ can be computed if}  \, n > 2k + (r-1) = %s .$'% latex(2*self.k+ dim_Sigma-1)))
#         if n_param > 2*self.k + (dim_Sigma - 1):
#             pretty_print(html(r'<p style="margin-top:2em; margin-bottom:2em; margin-left:4.5em">$$\mathbb{E}(%s) \; = \; %s $$</p>'  % (latex(new_E_inv_expr_lside) , latex(E_inv_num[Ik_indx-1]) )))

#     def evaluate_moment(self, Sigma, n_param, Ik_indx):
#         A = Sigma
#         dim_Sigma = Sigma.shape[0]

#         R2.<f,p,r> = QQ['f,p,r']

#         (lsideD, new_E_inv_expr_lside, new_E_inv_expr) = self.expectations_expressions(Ik_indx)

#         (Enum, E_inv_num)= self.wishart_expectations_numval(Sigma, n_param)

#         pretty_print(html(r'<div>$(i) = %s $</div>' % LatexExpr(self.P[Ik_indx-1])))
#         pretty_print(html(r'<p style= "margin-top:2em; margin-bottom:2em; margin-left:4.5em">$$2\Sigma = %s $$</p>' %latex(2*A) ))
#         pretty_print(html( r'<p style="margin-top:2em; margin-bottom:2em; margin-left:4.5em"> $$\mathbb{E}(%s) \; = \; %s$$</p>' % (latex(self.v_L[Ik_indx-1].subs(lsideD)/self.k) , latex(Enum[Ik_indx-1].subs({p:n/2})/self.k)) ))
#         pretty_print(html( r'$\text{The moments of } W^{-1} \text{ can be computed if}  \, n > 2k + (r-1) = %s .$'% latex(2*self.k+ dim_Sigma-1)))
#         if n_param > 2*self.k + (dim_Sigma - 1):
#             pretty_print(html(r'<p style="margin-top:2em; margin-bottom:2em; margin-left:4.5em">$$\mathbb{E}(%s) \; = \; %s $$</p>'  % (latex(new_E_inv_expr_lside) , latex(E_inv_num[Ik_indx-1]) )))

    def number_of_expectations(self):
        return self.number_of_partitions()

    def pretty_print_moment(self,t,inverse=False):
        m = self.moment(t,inverse)
        pretty_print(html(r'<p style="margin-top:2em; margin-bottom:2em; margin-left:4.5em">$ \mathbb{E}(%s) = %s $</p>' % (latex(m[0]),latex(m[1])) ))
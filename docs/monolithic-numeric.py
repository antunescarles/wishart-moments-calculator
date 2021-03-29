import numpy as np
from bisect import bisect_left ###

# Raising operator (R_ij^l(lambda))*
def raising(p,indexes,l):
    (i,j) = indexes
    mu = p
    mu[i]-=l
    mu[j]+=l
    mu[::-1].sort()
    return mu

# Busqueda del indice de la particion mu en la matriz de particiones.
def findIndex(mu, M, from_indx):
    for j in range(from_indx,M.shape[1]):
        if np.all(mu==M[:,j]) :
            return j

def jackCoefficients(k):
    
    verbose = False
    
    parts = list(Partitions(k).list())

    n = len(parts); #cardinal of the set of partitions

    ##### INITIALIZATION
    # Initialize the coefficient matrix with zeros in the valid positions and NaN in the invalid ones.
    C = np.zeros((n-1,n))
    
    #complete every partition with zeros and store them as columns in a matrix
    M = np.zeros((k,n))

    for i in range(0,n):
        partlen = len(parts[i])
        M[:partlen,i] = parts[i]

    if (verbose):
        print("M = \n", M, '\n')
        print("C = \n", C , '\n')

    ##### ALGORITHM

    #level of indentation (for formatting the output)
    level = 0

    #Calculation of the coefficients

    if (verbose):
        print("Calculo de los coeficientes \n")

    #At least for the moment, the partitions in the array are ordered in increasing lexicographic order, but the index of the greatest one is 0 and that of the last one is k-1. We may want to think how to change the index if we'd like.
    for s in range(0,n):
        level=1
        p = M[:,s]
        if (verbose):
            print(" "*level,"s = ", s," p = ", p)


        for i in range(0,k-1): #We need at least 2 indexes to do what we want.
            level=2
            if (verbose):
                print("  "*level,"i = ", i)


            for j in range(i+1,k):
                level=3
                if (verbose):
                    print("  "*level , "j = ", j)

                L= []

                if p[i]-p[j]>= 2 :
                    level=4

                    L = np.arange(1 , floor((p[i]-p[j])/2)+1)

                    if (verbose):
                        print("  "*level, "L = ", L)

                    for l in L:

                        level = 5

                        if (verbose):
                            print("  "*level, "l = ", l)
                        # raising() will change the argument inside its body. Array are mutable objects, so they are passed by reference in python functions. Be careful!
                        mu = np.array(p)
                        mu = raising(mu,(i,j),l)


                        #seach for the index of mu in the list of partitions
                        t = findIndex(mu, M, s)

                        if (verbose):
                            print("  "*level , "mu = ", mu, "--> t = ", t)

                        # assign the result (depending on the comparison between p[i]-l and p[j]+l)
                        # we could check wether the coefficient is already non-zero or not (if it is, we have already computed it and we can skip it).
                        if p[i]-l == p[j]+l:
                            level=6

                            # Count the number of ocurrences of p[i] - l 
                            n_pi_minus_l = np.count_nonzero(mu == (p[i]-l) )

                            #print("  "*level,"\n", C)

                            C[(n-1)-t,(n-1)-s] = (p[i]-p[j])*binomial(n_pi_minus_l , 2)
                            # The indexes are those because the partitions in the matrix M decrease as the index of the columns increase.
                            # i.e, mu[0] > mu[1] > ... > mu[n-1]
                            # This is more natural and efficient to compute the coefficients but it's not the way we want them to appear in the matrix to compute the Jack polynomials.

                            if (verbose):
                                print("  "*level, "caso igual ", "C[%d,%d]"  %(s,t), " = ",  C[(n-1)-t,(n-1)-s])
                        else:
                            level =6

                            # Count the number of ocurrences of p[i] - l and p[j] + l in mu
                            
                            n_pi_minus_l = np.count_nonzero(mu == (p[i]-l) )
                            n_pj_plus_l = np.count_nonzero(mu == (p[j]+l) )
                            C[(n-1)-t,(n-1)-s] = (p[i]-p[j])*n_pi_minus_l*n_pj_plus_l

                            if (verbose):
                                print("  "*level, "caso distinto ", "C[%d,%d]"  %(s,t), " = ",  C[(n-1)-t,(n-1)-s])
                else:
                    if (verbose):
                        print("  "*level, "L = ", L)
    return C

def jackDiagonal(k,s):
    
    partitions =Partitions(k).list() # This list is sorted in decreasing lexicographic order. So partitions[0] == [k].
    mu = partitions[s]   
    
    L = posets.IntegerPartitionsDominanceOrder(k);
    
    p = L.principal_order_ideal(mu) # principal_order_ideal returns a list of partitions/lists. 
                                               # Since lists in python are sorted by increasing lexicographic order by default,
                                               # we don't need to do more.
    
    R.<t>=ZZ['t']
    
    diag = []
    pp = [] # will become the list of completions of partitions p[i] (zero-padding)
    l_lam = [ len(p[i]) for i in range(0,len(p)) ]
    
    for i in range(0,len(p)):
        z = k-l_lam[i] # number of zeros we need to append to complete the partition
        pp.append(p[i]+z*[0])

    # Compute d_mu - d_lambda
    for i  in range(0,len(pp)-1): # the computatation stops in the len(pp)-1 bc we don't have to compute dmu-dmu. Remember pp[len(pp)-1] == mu <--> p[s]
        lam= pp[i]
        mu = pp[-1] # for a list l, l[-1] returns the last element of the list.
        a = 1
        if  l_lam[i] < l_lam[-1] :
            a = -1
            print("l_mu = %d"% l_mu, " l_lam[%d] = %d"%(i,l_lam[i]))
            
        diag.append(a*sum([(t/2)*(lam[j]^2 - mu[j]^2) - (j+1)*(lam[j]- mu [j]) for j in range(0,k) ]))
    return diag

def subdiagDet(J,k,subdiag,verbose):
    ## Computation of the determinant (Lapointe, Lascoux, Morse, 1999, Prop. 4)
    n = J.ncols()
    assert ( n == J.nrows()  ), "distinto numero de columans y filas en J"

    c = [-1]*n

    c[n-1] = prod([subdiag[i] for i in range(0,len(subdiag))])

    # Recursion

    i=n-2

    while i >= 0:
        c[i] = QQ['t'](sum([J[i+1,j]*c[j] for j in range(i+1,n)]))

        # Casteo xq J[_,_] es un polinomio en las variables m con coeficientes en QQ[t], no un elemento de QQ[t] directamente.
        # el -1  multiplicando a J[i+1,i] es porque en la formula hay que dividir por -a[i+1,i]

        c[i] = c[i]// QQ['t']((-1)*J[i+1,i])

        i -=1

    cc = [(factorial(k)*c[i] //c[0]) for i in range(0,n)]

    return cc

def computeJack(k,s,verbose):
    # Compute the jack polynomial coefficients in monomial and power-sum basis evaluated in t=2 (t is alfa in [Lapointe,Lascoux,Morse,1999])
    # returns a dictionary with two keys 'm' and 'p', corresponding to monomial and power-sum basis coeffs respectively.
    # The coeffs are given respect to the increasing order of partitions.
    # if verbose = True intermediate computations are displayed.

    n = Partitions(k).cardinality()

    assert (s>=0 ) , "s < 0 (here s = %d)" % s
    assert(s<= n-1) , "s > #partitions-1 (here s = %d and #partitions = %d )" % (s,n)

    C = jackCoefficients(k)

    diag = jackDiagonal(k,s)

    P = Partitions(k).list()

    ###

    ### List of variable names

    ## Find the indices of the partitions that are in the chain.

    # We create a fresh list with the partitions to inverse it. It's crucial tha it is sorted in increasing order to use bisect.
    Q = Partitions(k).list()
    Q.reverse()

    # print("P[(n-1)-s] = ",P[s],"== Q[s] = ",Q[(n-1)-s])

    L = posets.IntegerPartitionsDominanceOrder(k).principal_order_ideal(Q[(n-1)-s])

    M = len(L)# Remind that the list of partitions may exclude some of the partition succeded by mu in the dominance order. 
              # So there may not be exactly s partitions in the list.

    indices = [] # indices of the select partitions in the list of partitions sorted by increasing order, Q.
    for p in L:
        ix = bisect_left(Q,p) # We use bisect here because it's probably faster than l.index(x) specially if k is large.
        assert (ix < n), "element not found by bisect"
        indices.append(ix)
    # print("indices: ",indices)

    mvars = ''
    for i in range(1,M+1):
        if i!= 1:
            mvars= mvars+','
        mvars = mvars+ 'm_'+str(i)

    ## Try to implement this, for giving the variables a latex_name for printing them nicely
    # mvars = [var('m_'+str(i), latex_name = "m_{"+latex(L[i-1])+"}") for i in range(1,M+1)]
    # print("** ", mvars)

    R = QQ['t'][mvars]
    if verbose:
        print("variables definition string: ", mvars)

    #R.<t> = QQ['t'][]
    #J = matrix(SR,A)

    m_aux = R.gens()
    if verbose:
            print("variables (objects) m_aux: ", m_aux)

    if s != n-1:
        # WARNING: we need to select the corresponding columns of the matrix according to the partitions in the list 
        #(which are not necessarily all the partitions from 0 to s).

    #     A = np.array(C[0:n-s-1,0:n-s]) # Este paso es necesario, sino hay un error de contiguidad en memoria...
    #                                    # Remember that slicing start:end means [start,end).

        # We extract the submatrix corresponding to the partitions selected.
        # Ver si esto se puede hacer más prolijo indexando en los dos argumentos.
        AA = np.array(C[:,indices])
        A = np.array(AA[indices[0:len(indices)-1],:])

        J = matrix(R,A)

        # Fill the upper diagonal with the corresponding elements depending on s.

        for i in range(J.nrows()):
            if J.nrows() != len(diag):
                if verbose:
                    print("#J = %d y #diag = %d" % (J.nrows(),len(diag)),'\n')
            J[i,i]=diag[i]

        #m_aux = list(var('m_%d'%(i+1)) for i in range(J.ncols()))
        #m = matrix(SR, m_aux)

        m = matrix(R, m_aux)

        J = m.stack(J)

        if verbose:
            print("k = %d , s = %d" % (k,s))
            print('\n', J, '\n')

        ## Computation of the determinant (Lapointe, Lascoux, Morse, 1999, Prop. 4)

        cc = subdiagDet(J,k,diag,verbose)

        jpol = sum([cc[i]*m_aux[i] for i in range(0,len(cc))])

        # Output
        if verbose:
            print("Jack polynomial with [k = %d] and [s = %d] corresponding to partition" % (k,s),P[s], "\n")
            for j in range(0,M):
                print("(",cc[j].factor(),")","*",L[j])

        # Change to power-sum basis

        m = SymmetricFunctions(QQ['t']).monomial()
        p = SymmetricFunctions(QQ['t']).power()

        f= sum([cc[i](2)*m[L[i]] for i in range(0,M) ])

        if verbose:
            print("\nChange to power-sum basis (with t=2)\n")
            print("--MONOMIAL BASIS--")
            print(f,"\n")
            print("--POWER-SUM BASIS--")
            print(p(f),"\n")

        jcoefs = {} # empty dictionary
        # coefficients respect to the monomial basis (t=2)
        jcoefs['m'] = [cc[i](2) for i in range(0,len(cc))]

        # coeffiecients respect to the power-sum basis (t=2)
        P.reverse()
        #jcoefs['p'] = [p(f).coefficient(P[i]) for i in range(0,len(cc))]+[0]*s
        jcoefs['p'] = [p(f).coefficient(P[i]) for i in range(0,n)]

    #     print("jcoefs['m'] : ", jcoefs['m'])
    #     print("jcoefs['p'] : ", jcoefs['p'])

    if s == n-1:
        jpol = factorial(k)*m_aux[0]

        # Output
        if verbose:
            print("Jack polynomial with [k = %d] and [s = %d] corresponding to partition" % (k,s),P[s], "\n")
            print("(",jpol.monomial_coefficient(m_aux[0]),")","*",L[0])

        # Change to power-sum basis

        m = SymmetricFunctions(QQ['t']).monomial()
        p = SymmetricFunctions(QQ['t']).power()

        f= factorial(k)*m[L[-1]]

        if verbose:
            print("\nChange to power-sum basis (with t=2)\n")
            print("--MONOMIAL BASIS--")
            print(f,"\n")
            print("--POWER-SUM BASIS--")
            print(p(f),"\n")

        jcoefs = {} # empty dictionary
        # coefficients respect to the monomial basis (t=2)
        jcoefs['m'] = [factorial(k)]

        # coeffiecients respect to the power-sum basis (t=2)
        P.reverse()
        jcoefs['p'] = [p(f).coefficient(P[i]) for i in range(0,n)]#+[0]*s

    #     print("jcoefs['m'] : ", jcoefs['m'])
    #     print("jcoefs['p'] : ", jcoefs['p'])
    
    return jcoefs

######################################## PORTRAITS.py ####################################

def toPortrait(t,k):
    #  t is a type, or equivalently, a partition
    t = list(t) # we have to ensure we work with a list and not a object of another data type.
    
    i = [0]*k
    set_t = set(t)
    for j in set_t:
        #  we want to represent to store st such that st[0]*1 + st[1]*2 + st[2]*3 + ... + st[k-1]*k = k
        # notice that index starts from zero but is the same. That's the reason why we add 1 to i below:
        i[j-1] = list(t).count(j)
    return i

def traceDecorator(l,varname):
    # l sera j+1, la potencia del argumento
    # p sera i[j], la potencia de la traza
    a = "\\mathrm{tr}\\,"
    if (l == 1):
        a = a+ varname
    else:
        a = a+varname+"^%d"%(l)
    
    return "("+a+")"

def traceDecoratorInv(l,varname):
    # l sera j+1, la potencia del argumento
    # p sera i[j], la potencia de la traza
    a = "\\mathrm{tr}\\,"+varname+"^{-%d}"%(l)
    
    return "("+a+")"

def compute_r(i,k):
    # i is a portrait

    w = var('w')
    
    # When we have b1 we want tr\sigma instead of tr(\sigma^1)
#     r_i = prod([var('b%d'%(j+1),latex_name="\\mathrm{tr}(\\sigma^%d)"%(j+1))^(i[j]) for j in range(0,k) ])
    r_i = prod([var('b%d'%(j+1),latex_name = traceDecorator(j+1,"\\sigma") )^(i[j]) for j in range(0,k) ])
    return r_i

def compute_L(i,k):
    w = var('w')
    
    r_i = compute_r(i,k)
    
    L_i = sum([expand( r_i*(j+1)*i[j]*w^(j+1)/var('b%d'%(j+1)) ) for j in range(0,k) ])
    # ^  por alguna razon si multiplicamos r[i] afuera de sum([...]) no simplifica la bien la expresión...

    return L_i

def computeNum_r(i,k,S):
    tr = [ (S^(j+1)).trace() for j in range(0,k)]
    r_i = prod([ (tr[j])^(i[j]) for j in range(0,k) ])
    return (r_i , tr)

def computeNum_L(i,k,S):
    
    (r_i,tr) = computeNum_r(i,k,S)
    
    L_i = sum([ r_i*(j+1)*i[j]*S^(j+1)/tr[j] for j in range(0,k) ])
    # ^  por alguna razon si multiplicamos r[i] afuera de sum([...]) no simplifica la bien la expresión...

    return L_i

######################################### MAIN ######################################### 

r = A.nrows()
assert (r == A.ncols()), "The matrix given is not squared."
assert (A.is_symmetric()), "The matrix given is not symmetric."

outmost_verbose = False

## Check that the dim of the matrix is in fact a number of partitions.
i=1
r_correct = False #flag that indicates when r is indeed the number of partitions of some k.
passed = False # flag that is set true when we get to some i with Partitions(i).cardinality()> k

k=NaN
Ik_indx = 0

while (not r_correct and not passed and i<=10) :
    card = Partitions(i).cardinality()
#     print("r = %d and card = %d"% (r,card))
    if r == card:
        r_correct = True
        k = i
    else:
        if card > r:
            passed = True
    i+=1

assert (r_correct), "The size of the matrix given is not the cardinality of the number of partitions of any k (at least for k<=10)"

# TO CHANGE. Use r or n uniformly troughout the script.
n = r

s = 0 # partitions go like mu(n-(n-1)) = m(1) < mu(n-(n-2)) = mu(2) < mu(n-1) < mu(n-0) = [n]
      # So s can range from 0 to n-1
      # The program will compute the Jack polynomial corresponding to partition mu[n-s] of the list mu of partitions

print("k = %d , r = %d " % (k,n),"\n")



verbose = False # if verbose == True intermediate computations and checkings are displayed.



# Compute the jack polynomial of order k associated to partition s in the monomial and power basis, EVALUATED IN t=2
# coef = []
# coef.append(computeJack(k,0,False)['p'])
# coef.append(computeJack(k,1,False)['p'])
# coef.append(computeJack(k,2,False)['p'])
# coef.append(computeJack(k,3,False)['p'])
# coef.append(computeJack(k,4,False)['p'])
# coef.append(computeJack(k,5,False)['p'])
# coef.append(computeJack(k,6,False)['p']) # implementar caso s = n-1

coef = computeJack(k,s,verbose)
Bk = matrix(QQ,1,coef['p'])

t=0
#     print("s = %d " % 0,coef['p'])
t+=1
# for t in range(1,n-1):
while t<=n-1: # we use while instead of for bc when k=2 range(1,1) is empty and it never enters the loop
    coef = computeJack(k,t,verbose)
    row =  matrix(QQ,1,coef['p'])
#         print("s = %d " % t,coef['p'])
    Bk = Bk.stack(row)
    t+=1

if outmost_verbose: 
    print("\nBk : \n")
    print(Bk, "\n")

R2.<f,p> = QQ['f,p']

P = Partitions(k).list()

Dk = matrix(R2,n,n,0)



if outmost_verbose:  print("Elementos de la diagonal de Dk factorizados\n")
pm = [1]*n
for i in range(0,n):
    lm = len(P[i])
    for j in range(1,lm+1):
            for s in range(1,P[i][j-1]+1):
                pm[i] *= p +s-1- (j-1)*f
    Dk[i,i] = pm[i].subs({f:1/2}) # Evaluated in f = 1/2

    if outmost_verbose: print(P[i]," -->  ", pm[i].subs({f:1/2}).factor())

if outmost_verbose: 
    print("\nDk:\n")
    print(Dk) # Como mostrar las entradas de la matriz factorizadas? Aparentemente no hay una forma de hacerlo.

# Compute Mp 
IBk = Bk.inverse()

if outmost_verbose: 
    print("\n")
    print("IBk : \n")
    print(IBk,"\n")

Mp = IBk*Dk*Bk

if outmost_verbose: 
    print("Mp : \n")
    print(Mp)


## Computations of the moments

P.reverse()
# print(P,"\n")

r = []
L = []
Lnum = []

#     A = random_matrix(QQ,n,n)
#     A = identity_matrix(k)

for j in range(0,n):
    r.append(toPortrait(P[j],k))
    L.append(compute_L(r[j],k))
    Lnum.append(computeNum_L(r[j],k,2*A))

# for j in range(0,len(r)):
#     print(r[j]," ", compute_r(r[j]), " ", compute_L(r[j]))
# print("\nLnum: \n",Lnum)

# print("\nL: \n",L)

v_L = vector(SR,L)
# print(Lnum[0].parent())


#     print("\nv_Lnum: \n",v_Lnum)

if outmost_verbose: 
    print("\n")
    print("v_L = " , v_L,"\n")

    print("E = Mp*[r_(i)(w)]_(i): \n")

E = Mp*v_L
#     Enum = Mp*v_Lnum
Enum = [NaN]*n
#Para Enum hay que hacer las cuentas mas a mano porque no podemos formar un vector de matrices...
for i in range(0,n):
#     print(Mp[0,0]*Lnum[0])
    Enum[i] = sum([Mp[i,j]*Lnum[i] for j in range(0,n)])
#     print("\n",Enum[i],"\n")

if outmost_verbose: 
    for i in range(0,len(E)):
        print(E[i])

# Computations on demand
W = var('W')
N = var('N',latex_name="n")
S = var('S',latex_name="\\Sigma")

if outmost_verbose: print('\nPara Jero: \n')

#     Ik_indx = n-1
#     Ik_indx = n # Here we use indices starting at 1
if outmost_verbose: print("E[" ,v_L[Ik_indx-1].subs({w : W})/k ,"] = \n")

D = {p:N/2,w:2*S}

if outmost_verbose: print(E[Ik_indx-1].subs(D)/k,"\n")
show("2\\Sigma = "+ latex(2*A))
show("\\mathbb{E}("+latex(v_L[Ik_indx-1].subs({w : W})/k)+") = "+latex(Enum[Ik_indx-1].subs({p:N/2})/k))
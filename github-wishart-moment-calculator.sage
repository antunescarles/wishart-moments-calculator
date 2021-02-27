import numpy as np

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
    """

    for m in range(0,C.shape[1]):
        C[m:,m] = NaN
    """
    
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
    p=Partitions(k).list()
    R.<t>=ZZ['t']


    n = len(p)
    
    l=k-len(p[s])
    mu=[]
    cn=1 ## no se usa...

    diag = []
    

    for i in range(0,len(p[s])):
         mu.append(p[s][i])
    if l > 0:
        for i in range(1,l+1):
            mu.append(0)

    for i in range(1,len(p)-s):
        #print(i)
        r=len(p[len(p)-i])
        l=k-r
        lam=[]
        for j in range(0,r):
            lam.append(p[len(p)-i][j])
        if l > 0:
            for i2 in range(1,l+1):
                lam.append(0)  
        suma=0
        for j in range(0,k):
            suma=suma+(t/2)*(lam[j]^2-mu[j]^2)-(j+1)*(lam[j]-mu[j])
        ## 
        #print("suma = ", suma)
        diag.append(suma)

        #C[i-1,i-1] = suma
        #print("M[%d,%d] = " %(i-1,i-1), C[i-1,i-1], '\n')
    return diag

def computeJack(k,s,verbose):
    # Compute the jack polynomial coefficients in monomial and power-sum basis evaluated in t=2 (t is alfa in [Lapointe,Lascoux,Morse,1999])
    # returns a dictionary with two keys 'm' and 'p', corresponding to monomial and power-sum basis coeffs respectively.
    # The coeffs are given respect to the increasing order of partitions.
    # if verbose = True intermediate computations are displayed.
    
    n = Partitions(k).cardinality()

    assert (s>=0 ) , "s < 0 (here s = %d)" % s
    assert(s<= n-1) , "s > #partitions-1 (here s = %d and #partitions = %d )" % (s,n)
    #assert (s< n-1) , "case s = #partitions-1 not implemented (here s = %d and #partitions = %d )" % (s,n)


    C = jackCoefficients(k)

    diag = jackDiagonal(k,s)
    
    P = Partitions(k).list()

    ###

    # List of variable names
    mvars = ''
    #for i in range(1,n+1): ######################
    for i in range(1,n-s+1):
        if i!= 1:
            mvars= mvars+','
        mvars = mvars+ 'm_'+str(i)

    R = QQ['t'][mvars]
    if verbose:
        print("variables: ", mvars)

    #R.<t> = QQ['t'][]
    #J = matrix(SR,A)

    m_aux = R.gens()
    if verbose:
            print("m_aux: ", m_aux)

    if s != n-1:
        A = np.array(C[0:n-s-1,0:n-s]) # Este paso es necesario, sino hay un error de contiguidad en memoria...
                                       # Remember that slicing start:end means [start,end).
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

        #c = [-1]*n
        c = [-1]*(n-s)

        c[n-s-1] = -1*diag[0] # Remind in python indices start from 0.
        for i in range(1,len(diag)):
            c[n-s-1] = c[n-s-1]*(-1)*diag[i]

        # Recursion

        #i=n-2
        i = n-s-2

        while i >= 0:
            c[i] = QQ['t'](sum([J[i+1,j]*c[j] for j in range(i+1,n-s)]))

            # Casteo xq J[_,_] es un polinomio en las variables m con coeficientes en QQ[t], no un elemento de QQ[t] directamente.
            # el -1  multiplicando a J[i+1,i] es porque en la formula hay que dividir por -a[i+1,i]
            c[i] = c[i]// QQ['t']((-1)*J[i+1,i])

            i -=1

        cc = [(factorial(k)*c[i] //c[0]) for i in range(0,n-s)]

        jpol = sum([cc[i]*m_aux[i] for i in range(0,n-s)])

        

        # Output
        if verbose:
            print("Jack polynomial with [k = %d] and [s = %d] corresponding to partition" % (k,s),P[s], "\n")
            j=n-s
            while j > 0:
                print("(",cc[j-1].factor(),")","*",P[n-j])
                j -=1

        # Change to power-sum basis

        m = SymmetricFunctions(QQ['t']).monomial()
        p = SymmetricFunctions(QQ['t']).power()

        f= sum([cc[i](2)*m[P[n-1-i]] for i in range(0,n-s) ])

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

    if s == n-1:
        jpol = factorial(k)*m_aux[0]
        
        # Output
        if verbose:
            print("Jack polynomial with [k = %d] and [s = %d] corresponding to partition" % (k,s),P[s], "\n")
            print("(",jpol.monomial_coefficient(m_aux[0]),")","*",P[n-1])
            
        # Change to power-sum basis

        m = SymmetricFunctions(QQ['t']).monomial()
        p = SymmetricFunctions(QQ['t']).power()

        f= factorial(k)*m[P[n-1]]

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
        
    return jcoefs

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


## IMPORTANTE. Fijarse que cuando k=6 hay algo que falla al calcular los momentos.

#### PARAMETERS

# Given a partition mu, construct the matrix for the Jack polynomial associated to mu (related to the parameter s)
# It is basically, to extract a submatrix of C and add the diagonal elements to it

# k = 5
@interact
def _(k = input_box(3),s=input_box(0),Ik_indx = input_box(Partitions(3).cardinality())):

    outmost_verbose = False

    n = Partitions(k).cardinality()

    s = 0 # partitions go like mu(n-(n-1)) = m(1) < mu(n-(n-2)) = mu(2) < mu(n-1) < mu(n-0) = [n]
          # So s can range from 0 to n-1
          # The program will compute the Jack polynomial corresponding to partition mu[n-s] of the list mu of partitions

    print("k = %d , n = Partitions(k).cardinality() = %d , s = %d " % (k,n,s),"\n")

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
    for j in range(0,n):
        r.append(toPortrait(P[j],k))
        L.append(compute_L(r[j],k))

    # for j in range(0,len(r)):
    #     print(r[j]," ", compute_r(r[j]), " ", compute_L(r[j]))

    v_L = vector(SR,L)

    if outmost_verbose: 
        print("\n")
        print("v_L = " , v_L,"\n")

        print("E = Mp*[r_(i)(w)]_(i): \n")

    E = Mp*v_L
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

    for i in range(1,n+1):
#         D[var('b%d'%i)] = (2^i)*var('b%d'%i)
        D[var('b%d'%i)] = (2^i)*var('b%d'%i,latex_name = traceDecorator(i,"\\Sigma"))

    if outmost_verbose: print(E[Ik_indx-1].subs(D)/k,"\n")

    print("\n")

    show("\\mathbb{E}("+latex(v_L[Ik_indx-1].subs({w : W})/k)+") = "+latex(E[Ik_indx-1].subs(D)/k))
#     show("[\\mathbb{E}(L_{r_{(i)}}(U))]_{I_k} = "+latex(E[n-1].subs(D)/k))

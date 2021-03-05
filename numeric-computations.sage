import numpy as np

# I know we are advised against using this wildcards bc of the potential collisions of namespaces
# Just let's do it for now until we take the time to add the prefixes to the functions of this modules.
#from portraits import *
#from jacks import *
load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/jacks.sage')
load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/portraits.sage')

## IMPORTANTE. Fijarse que cuando k=6 hay algo que falla al calcular los momentos.

#### PARAMETERS

# Given a partition mu, construct the matrix for the Jack polynomial associated to mu (related to the parameter s)
# It is basically, to extract a submatrix of C and add the diagonal elements to it

# k = 5
@interact(layout = dict(top = [['k','s','Ik_indx']]) )
def _(k = input_box(3),s=input_box(0),Ik_indx = input_box(Partitions(3).cardinality())):

    outmost_verbose = True

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
    Lnum = []
    
#     A = random_matrix(QQ,n,n)
#     A = identity_matrix(k)
    
    for j in range(0,n):
        r.append(toPortrait(P[j],k))
        L.append(compute_L(r[j],k))
        Lnum.append(computeNum_L(r[j],k,2*A))

    # for j in range(0,len(r)):
    #     print(r[j]," ", compute_r(r[j]), " ", compute_L(r[j]))
    print("\nLnum: \n",Lnum)
    
    print("\nL: \n",L)

    v_L = vector(SR,L)
    print(Lnum[0].parent())
    
    
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
        print(Mp[0,0]*Lnum[0])
        Enum[i] = sum([Mp[i,j]*Lnum[i] for j in range(0,n)])
        print("\n",Enum[i],"\n")

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

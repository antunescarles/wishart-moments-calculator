import numpy as np

load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/jacks.sage')
load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/portraits.sage')

@interact
def wrpr(k = input_box(3,width = 8, label="$k$")):

    outmost_verbose = False

    assert (k >= 1) , "Error: k < 0"

    n = Partitions(k).cardinality()

    @interact
    # def _(k = input_box(3,width = 8, label="$k$"),Ik_indx = input_box(Partitions(3).cardinality(),width = 8, label ="$i$")):
    def _(Ik_indx = slider(1,n,step_size=1, label ="$i$")):

        s = 0 # partitions go like mu(n-(n-1)) = m(1) < mu(n-(n-2)) = mu(2) < mu(n-1) < mu(n-0) = [n]
        # So s can range from 0 to n-1
        # The program will compute the Jack polynomial corresponding to partition mu[n-s] of the list mu of partitions
        
        # Validation of the input
        assert (1 <= Ik_indx and Ik_indx <= n) , "Error: i < 0 or i > n (#partitions)"

        print("k = %d , r = Partitions(k).cardinality() = %d , s = %d " % (k,n,s),"\n")

        verbose = False # if verbose == True intermediate computations and checkings are displayed.

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

        ## Calculamos la otra diagonal para el momento de la inversa
        Dk_star = matrix(R2,n,n,0)

        # To do: Chequear condición....

        if outmost_verbose:  print("Elementos de la diagonal de Dk factorizados\n")
        R3.<q> =  QQ['q'];

        qm = [1]*n
        for i in range(0,n):
            lm = len(P[i])
            for j in range(k-lm+1,k+1):
                    for s in range(1,P[i][k-j+1 -1]+1):
                        qm[i] *= p + (k-j+1)*f -s # here I'd like to use another var, e.g, q instead of the same p,
                                                  # but as Ill inmediatelly substitute it's not worth the effort thinking a better solution.
            Dk_star[i,i] = qm[i].subs({f:1/2 , p : (p - n*f)}) # Evaluated in f = 1/2 and q= p-nf (o como aparece en el paper q = p-rf)

#             print(P[i]," -->  ", qm[i].subs({f:1/2}).factor())

        # if outmost_verbose: 
#         print("\nDk_star:\n")
#         print(Dk_star)

        # When it corresponds, compute M^*(p-rf) r = Partitions(k).cardinality() == n
        M_pnf_star = IBk*Dk_star*Bk

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
        Sinv = var('Sinv', latex_name = "\\Sigma")
        Winv = var('Winv',latex_name = "W^{-1}")

    #     Ik_indx = n-1
    #     Ik_indx = n # Here we use indices starting at 1
        if outmost_verbose: print("E[" ,v_L[Ik_indx-1].subs({w : W})/k ,"] = \n")

        D = {p:N/2,w:2*S}
        Dinv = {p:N/2 , w : 2*Sinv^(-1) }

        for i in range(1,n+1):
    #         D[var('b%d'%i)] = (2^i)*var('b%d'%i)
            D[var('b%d'%i)] = (2^i)*var('b%d'%i,latex_name = traceDecorator(i,"\\Sigma"))

            Dinv[var('b%d'%i)] = (2^(-i))*var('a%d'%i,latex_name = traceDecoratorInv(i,"{\\Sigma}"))


        if outmost_verbose: print(E[Ik_indx-1].subs(D)/k,"\n")

#         print("\n")
#         print(E[Ik_indx-1].subs(D)/k,"\n")
        E_inv_expr = E[Ik_indx-1].subs(Dinv)/k
#         print(E_inv_expr,"\n")

        ## Artifact to print E[\Sigma ^{-1}] nicely
        # 1) Extract the coefficients of every negatice power of Sinv
        # 2) Form a new expression multiplying the coef of the (-j)-th powe of Sinv for a new variable, something like Sj with latex_name \Sigma^{-j}
        
        l = E_inv_expr.coefficients(Sinv)
#         print("list of coeff of Sinv with its exponents: \n")
#         print(l)
        new_E_inv_expr = sum( [ c[0]*var('S%d'%(-c[1]), latex_name = "{\\Sigma^{%d}}"%c[1]) for c in l] )
#         print("New expression: \n")
#         print(new_E_inv_expr)
#         show(latex(new_E_inv_expr))

        print("\n")

        lsideD = {w : W}
        lsideDinv = {w:Winv}

        for i in range(1,n+1):

            lsideD[var('b%d'%i)] = var('b%d'%i,latex_name = traceDecorator(i,"{W}"))

            lsideDinv[var('b%d'%i)] = var('a%d'%i,latex_name = traceDecoratorInv(i,"{W}"))

        show("\\mathbb{E}("+latex(v_L[Ik_indx-1].subs(lsideD)/k)+") \\; = \\; " +latex(E[Ik_indx-1].subs(D)/k))
        print("\n")
        # Show expectation of the inverses
        show("\\text{And if }\\, n \\geq 2k + (r-1) =  " + latex(2*k + n - 1))
        show("\\mathbb{E}("+latex(v_L[Ik_indx-1].subs(lsideDinv)/k)+") \\; = \\; "+latex(new_E_inv_expr))

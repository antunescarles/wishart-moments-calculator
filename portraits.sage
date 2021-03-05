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

def computeNum_r(i,k,S):
    tr = [ (S^(j+1)).trace() for j in range(0,k)]
    r_i = prod([ (tr[j])^(i[j]) for j in range(0,k) ])
    return (r_i , tr)

def computeNum_L(i,k,S):
    
    (r_i,tr) = computeNum_r(i,k,S)
    
    L_i = sum([ r_i*(j+1)*i[j]*S^(j+1)/tr[j] for j in range(0,k) ])
    # ^  por alguna razon si multiplicamos r[i] afuera de sum([...]) no simplifica la bien la expresión...

    return L_i
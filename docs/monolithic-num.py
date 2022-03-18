load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/sagecell/ObjectWithPartitions.py')
load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/sagecell/Jacks.py')
load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/sagecell/Expectations.py')

############### Interactive interface ###########################

@interact
def wrpr(k = input_box(2,width = 8, label="$k$") , N_param = input_box(2,width = 8,label = "$n$") , positive =checkbox(False,"Compute moment of $W$"), inverse = checkbox(True,label = 'Compute moment of $W^{-1}$') ):

    outmost_verbose = False

    assert (k >= 1) , "Error: k < 0"

    wishartk = Expectations(k)
    
    i_list = [1 .. wishartk.number_of_expectations()]
    i_list.reverse()

    @interact
    def _(Ik_indx = slider(vmin = i_list, label ="$(i)$")):
        
        # Validation of the input
        assert (1 <= Ik_indx and Ik_indx <= wishartk.number_of_expectations()) , "Error: i < 0 or i > n (#partitions)"
        
        pretty_print(html( r'$(i) = %s$' %  latex(tuple(wishartk.partition_to_portrait(wishartk.P[Ik_indx-1]))) ))

        if positive:
            wishartk.pretty_print_eval_moment(Ik_indx-1,N_param,Sigma)

        if inverse:
            if (N_param > 2*wishartk.k + Sigma.shape[0] -1):
                wishartk.pretty_print_eval_moment(Ik_indx-1,N_param,Sigma,True)
            else:
                pretty_print(html( r'$  \text{ The integer } n \text{ must satisfy } n > %d \text{ to compute }  \mathbb{E}(%s) $' % ((2*wishartk.k + Sigma.shape[0] -1),latex(wishartk.expression(Ik_indx-1,True)))))

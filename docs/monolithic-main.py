load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/sagecell/ObjectWithPartitions.py')
load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/sagecell/Jacks.py')
load('https://raw.githubusercontent.com/antunescarles/wishart-moments-calculator/main/sagecell/Expectations.py')


############### Interactive interface ##########################33


@interact
def wrpr(k = input_box(2,width = 8, label="$k$")):

    outmost_verbose = False

    assert (k >= 1) , "Error: k < 0"

    
    wishartk = Expectations(k)
    
    i_list = [1 .. wishartk.number_of_expectations()]

    i_list.reverse()

    @interact
    def _(Ik_indx = slider(vmin = i_list, label ="$(i)$"), positive =checkbox(False,"Compute moment of $W$"), inverse = checkbox(False,label = 'Compute moment of $W^{-1}$') ):

        # Validation of the input
        assert (1 <= Ik_indx and Ik_indx <= wishartk.number_of_expectations()) , "Error: i < 0 or i > n (#partitions)"
        
        pretty_print(html( r'$(i) = %s$' %  latex(tuple(wishartk.partition_to_portrait(wishartk.P[Ik_indx-1]))) ))

        if positive:
            wishartk.pretty_print_moment(Ik_indx-1)

        if inverse:
            pretty_print(html( r'$\text{If } \, n > %d + (r-1)$' % (2*wishartk.k)))
            wishartk.pretty_print_moment(Ik_indx-1,True)
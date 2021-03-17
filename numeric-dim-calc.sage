@interact
def _(k = input_box(3,width = 5, label="$k$")):
    
    assert (k >= 1) , "Error: k < 0"

    n = Partitions(k).cardinality()

    print("A matrix A of dimension %dx%d is expexted \n" % (n,n))

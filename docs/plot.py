import matplotlib.pyplot as plt

def plot(r,res,theo=None):
    _,ax1 = plt.subplots()
    l1 = ax1.plot(r,res,'bx',label="Code results")

    if theo is not None:
        ax2 = ax1.twinx()
        l2 = ax2.plot(r,theo,'r',label="Expected results")
        lns = l1+l2
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc=0)
    else:
        ax1.legend(loc=0)

    # added these three lines

    plt.grid()
    plt.show()
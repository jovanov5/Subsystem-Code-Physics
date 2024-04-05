# import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('Agg')

def edge_picker(i, j, direction, system):
    """
    :param i: int, the x coordinate of the vertex
    :param j: int, the y coordinate of the vertex
    :param direction: int, the direction of the edge, 0 for horizontal, 1 for vertical
    """
    return 2*(i+j*system.L)+direction

def visualise_the_stabiliser(stabiliser, system):
    stab = list(stabiliser) # This is given as a 2nbit long touple that needs to be make into a list
    L = system.L
    p = plt.figure().add_subplot(111)
    nbits = system.nbits

    for i in range(L):
        for j in range(L):
            if stab[edge_picker(i, j, 1, system)] == 1:
                p.plot([i, i+1], [j, j], color='blue')
                p.plot([i+0.5, i+0.5], [j+0.5, j-0.5], color='blue', linestyle='dashed')
            if stab[edge_picker(i, j, 1, system)+nbits] == 1:
                p.plot([i, i+1], [j, j], color='red')
            if stab[edge_picker(i, j, 0, system)] == 1:
                p.plot([i, i], [j, j+1], color='blue')
                p.plot([i+0.5, i-0.5], [j+0.5, j+0.5], color='blue', linestyle='dashed')
            if stab[edge_picker(i, j, 0, system)+nbits] == 1:
                p.plot([i, i], [j, j+1], color='red')

    p.set_xlim(-1, L+1)
    p.set_ylim(-1, L+1)
    plt.show()


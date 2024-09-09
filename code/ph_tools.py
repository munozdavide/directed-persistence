import numpy as np
import matplotlib.pyplot as plt
import dionysus as dio


def persistence_dgms(filtration, max_dim=1):
    '''
    This function computes persistence diagrams from Gudhi simplex tree.

    Args:
        filtration (gudhi.SimplexTree): simplicial filtration as a Gudhi simplex tree
        max_dim (int): maximum dimension to compute persistence; default is 1

    Returns:
        dgms (numpy.array): 2d array of persistence points ()
    '''

    dgms = filtration.persistence()

    dims = np.array([dgms[i][0] for i in range(len(dgms))])
    births, deaths = np.array([dgms[i][1] for i in range(len(dgms))]).T

    return np.array([dims, births, deaths])

def persistence_dgms_from_dionysus(filtration, max_dim=1):
    '''
    This function computes persistence diagrams from Dionysus filtration.

    Args:
        filtration (dionysus.Filtration): Dionysus filtration
        max_dim (int): maximum dimension to compute persistence; default is 1

    Returns:
        dgms (numpy.array): 2d array of persistence points
    '''

    dgms_dio = dio.init_diagrams(dio.homology_persistence(filtration), filtration)

    dims, births, deaths = [], [], []
    for i,dgm in enumerate(dgms_dio[:max_dim+1]):
        dims = dims+[i for _ in dgm]
        births = births+[pt.birth for pt in dgm]
        deaths = deaths+[pt.death for pt in dgm]
    births = np.array(births)
    deaths = np.array(deaths)

    return np.array([dims, births, deaths])

def plot_dgms(dgms, figtitle=None, filename=None, print_repeats=False,
              max_dim=1, ax=None, inf_val=0, alpha=0.3):
    '''
    Plots persistence diagrams
    '''

    # Colors for dimensions, up to 4
    colors = ['red', 'blue', 'green', 'orange', 'brown']

    from copy import copy
    dims, births, deaths = copy(dgms)

    max_death = max(deaths[deaths<np.inf])
    max_birth = max(births[births<np.inf])
    inf_val = max(1.05*max(max_death, max_birth), inf_val)

    deaths[ deaths==np.inf ] = inf_val

    if ax == None: ax = plt.gca()

    if figtitle != None: ax.set_title(figtitle)

    # Set proportional axes limits
    ax.set_xlim((-0.03*inf_val,inf_val))
    ax.set_ylim((-0.03*inf_val,inf_val*1.03))

    # Plot diagonal
    ax.plot([-0.05*inf_val,inf_val],[-0.05*inf_val,inf_val], 'k--', alpha=0.4)

    # Plot infinity line
    ax.plot([-0.1*inf_val,1.1*inf_val],[inf_val,inf_val], 'k-', alpha=0.3, label='\u221e')

    ax.grid(alpha = 0.2)

    # Plot every point
    visible_dims = [False]*5
    for i in range(len(dims)):
        ax.plot([births[i],], [deaths[i],], 'o', c=colors[int(dims[i])], alpha=alpha)
        visible_dims[int(dims[i])] = True

    # Plot hidden points to add labels
    for i,d in enumerate(visible_dims):
        if d: ax.plot([-1,],[-1,], 'o', c=colors[i], label='Dim '+str(i))
    ax.legend(loc='lower right')

    if filename != None: ax.savefig(filename+'.png', dpi=400)

    if print_repeats:
        count_dict = {i:dgms.count(i) for i in dgms}
        counts = np.array(list(count_dict.values()))
        repeats = np.array(list(count_dict.keys()))[counts > 1]
        print('---')
        for r in repeats:
            print('Point', tuple(r), 'shows up', count_dict[tuple(r)], 'times')

    return ax
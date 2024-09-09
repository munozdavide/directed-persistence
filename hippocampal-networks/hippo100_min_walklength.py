import sys
sys.path.append('../code')

from ph_tools import *
from walklength import *



for numHoles in range(5):
    numTrial = 101

    while numTrial <= 120:
        import scipy.io
        mat = scipy.io.loadmat('../../simHippo/allDissim/t_dissim_'+str(numHoles)+'_'+str(numTrial)[1:]+'.mat')
        mat = mat['dissim']

        # Shift edge weights to start at zero
        mat -= np.min(mat)

        print('Working on network', numHoles, str(numTrial)[1:], 'with size', mat.shape[0], end='\r', flush=True)

        f = walklength_filtration(mat.T)
        dgms = persistence_dgms(f)

        plt.figure(figsize=(5,5))
        plot_dgms(dgms, alpha=0.3, inf_val=4.2);
        plt.savefig('../../figures/hippocampal-networks-diagrams/walklen_min_dgm_'+str(numHoles)+'_'+str(numTrial)[1:]+'_'+str(mat.shape[0])+'.png', dpi=400);
        plt.close()

        dgms = dgms.T

        filename = '../../hippocampal-networks/dowker/walklen_min_dgms_hipponet_'+str(numHoles)+'_'+str(numTrial)[1:]+'.csv'

        import csv
        with open(filename, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(dgms_as_np)
        file.close()

        numTrial += 1


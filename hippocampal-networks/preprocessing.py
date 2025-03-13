# Function for selecting preprocessing for simulated hippocampal data

# Options are
# preprocs = ['original','min1', 'minpurge', 'purge','originalshort','min1short', 'minpurgeshort', 'purgeshort']

import numpy as np
from scipy.sparse.csgraph import shortest_path

def process_weights(weight_mat, method='original'):
	'''
	This function modifies the weights in a given network by the given preprocessing method.

	Args:
		weight_mat (numpy.array): Square 2d array where entry [i,j] is the weight of directed edge (i,j).
		method (str): Preprocessing method in {'original','min1', 'minpurge', 'purge','originalshort','min1short', 'minpurgeshort', 'purgeshort'}

	Returns:
		weight_mat (numpy.array): Array with modified weights
	'''

	if method not in ['original','min1', 'minpurge', 'purge','originalshort','min1short', 'minpurgeshort', 'purgeshort']:
		print('Invalid method')
		return None

	if 'purge' in method:
		# Purge unwanted edges
		weight_mat[weight_mat==np.max(weight_mat)] = np.inf

	if 'min' in method:
		# Shift edge weights to start at zero
		weight_mat -= np.min(weight_mat)
		if 'min1' in method:
			# Put max weights back on 1
			weight_mat[weight_mat==np.max(weight_mat)] = 1

	if 'short' in method:
		# Use shortest-path network
		weight_mat = shortest_path(weight_mat.T)

	return weight_mat
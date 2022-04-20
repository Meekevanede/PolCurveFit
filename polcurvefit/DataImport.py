import numpy as np

# Function to load an example file

def read_nova(filename):
	"""
		Function to import data from a nova file (Metrohm autolab)

		:param filename: path to the data file
		:type filename: string

		:rparam filename: 2 N-length arrays, containing the polarization curve data (E, I)
		:rtype filename: tuple
	"""
	
    data = np.loadtxt(filename,comments = '#', delimiter = '\t')
    return data[:,0], data[:,1] 

def read_gamry(filename):
	"""
		Function to import data from a gamry file (Gamry potentiostat)

		:param filename: path to the data file
		:type filename: string

		:rparam filename: 2 N-length arrays, containing the polarization curve data (E, I)
		:rtype filename: tuple
	"""

	data = np.loadtxt(filename,comments = '#', delimiter = '\t')
    return data[:,0], data[:,1]


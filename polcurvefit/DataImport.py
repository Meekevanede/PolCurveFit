import numpy as np

# Functions to load potential curve data

def load_txt(filename, lines_header = 0, columns_E_I = [0,1]):
	"""
		Function to import data from a generic text file, columns separated by space or tab.

		:param filename: path containing the data file
		:type filename: string

		:param lines_header: The number of lines of header information (Default:0) 
		:type lines_header: int

		:param columns_E_I: The columns containing the potential data E and the current data I. (Default: the first column contains E, second column I,[0,1])
		:type columns_E_I: 2-length sequence

		:returns: 2 N-length arrays, containing the polarization curve data (E, I)
		:rtype: tuple
	"""
	data = np.loadtxt(filename, skiprows = lines_header, usecols= (columns_E_I[0],columns_E_I[1]))
	return data[:,0], data[:,1] 

def load_nova(filename):
	"""
		Function to import data from a text file exported by the nova software (Metrohm autolab)

		:param filename: path containing the data file
		:type filename: string

		:returns: 2 N-length arrays, containing the polarization curve data (E, I)
		:rtype: tuple
	"""
	data = np.loadtxt(filename, skiprows = 1, delimiter = ';', usecols = (0,2))
	return data[:,0], data[:,1] 

def load_gamry(filename):
	"""
		Function to import data from a text (DTA format) file exported by gamry software (Gamry)

		:param filename: path containing the data file
		:type filename: string

		:returns: 2 N-length arrays, containing the polarization curve data (E, I)
		:rtype: tuple
	"""

	file = open(filename)
	counter = 1
	for line in file:
		if line.startswith("\t#"):
			break
		counter += 1

	data = np.loadtxt(filename, skiprows= counter, delimiter = '\t', usecols = (3,4))
	return data[:,0], data[:,1]


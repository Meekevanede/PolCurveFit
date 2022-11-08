import os
import numpy as np
import pandas as pd
import math
from scipy.optimize import curve_fit
from polcurvefit.forward import *
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib as mpl
import matplotlib.colors as colours

class polcurvefit:
	
	"""
	polcurvefit is a python class that can be used to analyze measured polarization curves and obtain parameters such as the corrosion potential, 
	Tafel slopes, corrosion current density and exchange current densities. The data can be fitted with 3 techniques: Tafel extrapolation (linear fit), 
	'Activation control fit' & 'mixed activation-diffusion control fit'. 

	:param E: The electrical potentials [V vs ref]
	:type E: N-length sequence

	:param I: The electrical current or current density [A or A/area]
	:type I: N-length sequence

	:param R: The resistance, to correct for the voltage drop (IR-drop/ohmic potential drop) [Ohm] (Default = 0.0)
	:type R: float

	:param sample_surface: The surface area of the metal sample surface measured. This used to convert from current to current density. If I is already the current density, use the default. If the surface area is not given, the obtained corrosion rates and exchange current density are simply currents (Default = 1.0).
	:type sample_surface: float

	"""

	def __init__(self,E,I, R = 0.0, sample_surface = 1.0):

		# check shape
		if len(np.array(E).shape) != 1:
			raise ValueError("The input potential should be a 1 dimensional sequence")
		if len(np.array(I).shape) != 1:
			raise ValueError("The input current should be a 1 dimensional sequence")

		# if input is a list, check if it contains only numbers
		if all(isinstance(e, (int, float)) for e in E) == False:
			raise ValueError("non-real numbers in potential input")
		if all(isinstance(e, (int, float)) for e in I) == False:
			raise ValueError("non-real numbers in current input")

		# Delete any zero values of I
		zero_indices = np.where(I == 0.0)[0]
		I = np.delete(I, zero_indices)
		E = np.delete(E, zero_indices)

		self.E = E
		self.I = I

		if E[0]>E[1]:
			self.E = np.flip(E)
			self.I = np.flip(I)

		self.E_obs = self.E

		# IR correction
		self.E = self._ir_correction(R)

		# obtain current density
		self.i = self.I/sample_surface

		# Initializing fitting parameters and results
		self.fit_results = None
		self.E_corr = None
		self.I_corr = None  	# This is the 10 base logarithmic of I_corr
		self.b_a = None
		self.b_c = None
		self.i_L = None         # This is the 10 base logarithmic of I_L
		self.gamma = None
		self.RMSE = None
		self.io_an = None
		self.io_cath = None
		self.w_ac = None
		self.W = None
		self.apply_W = None

		# Initializing parameters and results of the sensitivity analysis
		self.w_dc = None
		self.window_sens = None
		self.df_sens = None
		self.df_sens_std = None

	def __str__(self):
		return str(self.__class__) + ": " + str(self.__dict__)

		
############################## Tafel extrapolation (Linear fit)  ###################
####################################################################################
	def linear_fit(self, window, E_corr = 0.0, obtain_io = False, E_rev = 0):
		"""
		Fitting of a linear line to the data. This option is suitable for data which is solely under activation control on the anodic 
		or on the cathodic branch. It fits a linear line the data in the selected window. If the corrosion potential E_corr is given,
		the returned I_corr relects the corrosion rate.
		
		:param window: Lower and upper bounds of the data to be fitted [min value, max value], units [V vs ref] 
		:type window: 2-length sequence

		:param E_corr: The corrosion potential (OCP) [V vs ref]. If specified, the returned I_corr_best reflects the corrosion rate. (Default: 0.0)
		:type E_corr: float

		:param obtain_io: If True, the (cathodic OR anodic) exchange current density will be determined from the obtained Tafel slope and E_rev. (Default: False)
		:type obtain_io: bool

		:param E_rev: The reversible potential [V vs ref] of the anodic OR cathodic reaction of the Tafel region that is fitted. Needs to be specified if obtain_io=True. (Default = 0.0)
		:type E_rev: float

		:return tuple Results:  
			- fit_results (:py:class:`2xN array`) - the fit to the data [current densities (N-array), potentials (N-array)]
			- E_corr (:py:class:`float`) - the corrosion potential [V vs ref]
			- I_corr (:py:class:`float`) - the corrosion rate (if the corrosion potential is given as input) [A/surface area]
			- b (:py:class:`float`) - the obtained Tafel slope [V]
			- RMSE (:py:class:`float`) - the root mean squared error of the fit
			- io (:py:class:`float`) - Only returned, if obtain_io is True. The exchange current density [A/surface area].
		"""

		# select data in specified window
		E_cut, I_cut = self._windowcut(self.E,self.i,0.0, window)

		# fitting
		xdata = np.zeros((len(E_cut),1))
		xdata[:,0] = E_cut
		popt, pcov = curve_fit(linear_fit, E_cut, np.log10(abs(np.array(I_cut))))
		d_final = 10**np.array(linear_fit(E_cut, popt[0], popt[1]))
		if popt[0]<0:
			d_final = -1*d_final
		
		# Compute root mean squared error
		RMSE=math.sqrt(np.sum(np.square(I_cut-d_final)).mean())

		# Obtain real corrosion rate (if E_corr is given as input)
		I_corr_best = popt[0] * E_corr + popt[1] 
		
		self.fit_results = [d_final, E_cut]
		self.E_corr = E_corr
		self.I_corr = I_corr_best
		self.b_a = 1/popt[0]
		self.b_c = None
		self.RMSE = RMSE
		self.i_L = None

		# compute exchange current density	
		if obtain_io:
			io = self._compute_exchcurrent(I_corr_best,E_corr, 1/popt[0], E_rev)
			self.io_an = io
			return [d_final, E_cut], E_corr, 10**I_corr_best, 1/popt[0], RMSE, io
		else:
			return [d_final, E_cut], E_corr, 10**I_corr_best, 1/popt[0], RMSE


############################## Activation control fit ########################
##############################################################################
	def active_pol_fit(self,window, i_corr_guess = None, obtain_io = False, E_rev_an = 0, E_rev_cath = 0):
		"""
		Fitting of a theoretical description, representative for the anodic and cathodic activation controlled currents around the Corrosion potential (OCP).
		This option is suitable for data of which the anodic and cathodic branches are solely under activation control.
		
		:param window: Lower and upper bounds of the data to be fitted, relative to the corrosion potential [min value, max value], units [V vs E_corr] 
		:type window: 2-length sequence

		:param i_corr_guess: First guess of (the order of magnitude of) the corrosion current density (optional) [A/surface area], it might lead to faster convergence (Default: 10^-2).
		:type i_corr_guess: float

		:param obtain_io: If True, the exchange current densities will be determined from the obtained Tafel slope and E_rev_an & E_rev_cath. (Default: False)
		:type obtain_io: bool 

		:param E_rev_an: The reversible potential [V vs ref] of the anodic reaction of the data that is fitted. Needs to be specified if obtain_io=True. (Default = 0.0)
		:type E_rev_an: float

		:param E_rev_cath: The reversible potential [V vs ref] of the cathodic reaction of the data that is fitted. Needs to be specified if obtain_io=True. (Default = 0.0)
		:type E_rev_cath: float
		
		:return tuple Results: 
			- fit_results (:py:class:`2xN array`) - the fit to the data [current densities (N-array), potentials (N-array)]
			- E_corr (:py:class:`float`) - the corrosion potential [V vs ref]
			- I_corr (:py:class:`float`) - the corrosion rate  [A/surface area]
			- b_an (:py:class:`float`) - the obtained anodic Tafel slope [V]
			- b_cath (:py:class:`float`) - the obtained cathodic Tafel slope [V]
			- RMSE (:py:class:`float`) - the root mean squared error of the fit
			- io_an (:py:class:`float`) - Only returned, if obtain_io is True. The anodic exchange current density [A/surface area].
			- io_cath (:py:class:`float`) - Only returned, if obtain_io is True. The cathodic exchange current density [A/surface area].
			 
		"""
		if i_corr_guess == None:
			i_corr_guess = np.mean(self.i)

		# check input 
		if (i_corr_guess>self.i.max() or i_corr_guess<self.i.min()):
			raise ValueError("Specified i_corr_guess does not lie within the range of the current density")

		# obtain E_corr
		E_corr = self._find_Ecorr()

		# select data in specified window
		E_cut, I_cut = self._windowcut(self.E,self.i,E_corr, window)

		# fitting
		xdata = np.zeros((len(E_cut),2))
		xdata[:,0] = E_cut
		xdata[:,1] = E_corr

		popt, pcov = curve_fit(forward_normal, xdata, I_cut, p0=[math.log10(np.abs(i_corr_guess)), 0.0600, 0.100], bounds=((math.log10(np.abs(self.i).min()),0,0),(math.log10(np.abs(self.i).max()),1,1)))
		d_final = forward_normal(xdata, popt[0], popt[1], popt[2])
		d_final[d_final == 0] = 1E-8 # to get rid of zeros
		RMSE=math.sqrt(np.sum(np.square(I_cut-d_final)).mean())
		
		self.fit_results = [d_final, E_cut]
		self.E_corr = E_corr
		self.I_corr = popt[0]
		self.b_a = popt[1]
		self.b_c = popt[2]
		self.RMSE = RMSE
		self.i_L = None


		if obtain_io:
			io_an, io_cath = self._compute_exchcurrent(popt[0],E_corr, popt[1], E_rev_an, popt[2], E_rev_cath)
			self.io_an = io_an
			self.io_cath = io_cath
			return [d_final, E_cut], E_corr, 10**popt[0], popt[1], -popt[2], RMSE, io_an, io_cath
		else:
			self.io_an = None
			self.io_cath = None
			return [d_final, E_cut], E_corr, 10**popt[0], popt[1], -popt[2], RMSE


############################## Mixed activation-difussion control fit #############
###################################################################################
	def mixed_pol_fit(self,window, i_corr_guess = None, i_L_guess = None, fix_i_L = False, apply_weight_distribution = False, w_ac = 0.04, W = 75, obtain_io = False, E_rev_an = 0, E_rev_cath = 0):
		
		"""
		Fitting of a theoretical description, representative for the anodic activation controlled currents and cathodic mixed activation-diffusion controlled currents.
		This option is suitable for data of which the anodic is under activation control and the cathodic branche are under mixed activation-diffusion control.
		To obtain more accurate and less subjective fit, a specific weight distribution can be applied (see Documentation)
		
		:param window: Lower and upper bounds of the data to be fitted, relative to the corrosion potential [min value, max value], units [V vs E_corr] 
		:type window: 2-length sequence

		:param i_corr_guess: First guess of (the order of magnitude of) the corrosion current density (optional) [A/surface area], it might lead to faster convergence (Default: 10^-2).
		:type i_corr_guess: float
		
		:param i_L_guess: First guess of the limiting current density (optional) [A/surface area], it can to faster convergence and a more accurate fit (Default: 10^1).
		:type i_L_guess: float
		
		:param fix_i_L: If True, the limiting current density is fixed at the value of i_L_guess. (Default:False)
		:type fix_i_L: bool

		:param apply_weight_distribution: If True, the data within w_ac are given a certain percentage, W, of the overall weight. The application of this weight distribution leads generally to more consistent and accurate results.
		:type apply_weight_distribution: bool

		:param w_ac: determines the window activation control: the window around corrosion potential (E_corr) [-w_ac,+w_ac [V vs E_corr]], in which the data are given a certain weight percentage, W, of the overall weight of the data.
		:type w_ac: float
		
		:param W: The percentage of the weight assigned to the data within +/-w_ac around the corrosion potential. 
		:type W: float

		:param obtain_io: If True, the exchange current densities will be determined from the obtained Tafel slope and E_rev_an & E_rev_cath. (Default: False)
		:type obtain_io: bool 

		:param E_rev_an: The reversible potential [V vs ref] of the anodic reaction of the data that is fitted. Needs to be specified if obtain_io=True. (Default = 0.0)
		:type E_rev_an: float

		:param E_rev_cath: The reversible potential [V vs ref] of the cathodic reaction of the data that is fitted. Needs to be specified if obtain_io=True. (Default = 0.0)
		:type E_rev_cath: float

		:return tulple Results: 
			- fit_results (:py:class:`2xN array`) - the fit to the data [current densities (N-array), potentials (N-array)]
			- E_corr (:py:class:`float`) - the corrosion potential [V vs ref]
			- I_corr (:py:class:`float`) - the corrosion rate  [A/surface area]
			- b_an (:py:class:`float`) - the obtained anodic Tafel slope [V]
			- b_cath (:py:class:`float`) - the obtained cathodic Tafel slope [V]
			- i_L (:py:class:`float`) - the obtained limiting current density [A/surface area]
			- RMSE (:py:class:`float`) - the root mean squared error of the fit
			- gamma (:py:class:`int`) - Curvature defining constant (see Documentation)
			- io_an (:py:class:`float`) - Only returned, if obtain_io is True. The anodic exchange current density [A/surface area].
			- io_cath (:py:class:`float`) - Only returned, if obtain_io is True. The cathodic exchange current density [A/surface area].

		"""

		# obtaining automatic guesses for i_icorr and i_L
		if i_corr_guess == None:
			i_corr_guess = np.mean(self.i)
		
		if i_L_guess == None:
			i_L_guess = np.mean(self.i)

		# check input 
		if (i_corr_guess>self.i.max() or i_corr_guess<self.i.min()):
			raise ValueError("Specified i_corr_guess does not lie within the range of the current density")
		if (i_L_guess>self.i.max() or i_L_guess<self.i.min()):
			raise ValueError("Specified i_L_guess does not lie within the range of the current density")

		# obtain E_corr
		E_corr = self._find_Ecorr()

		# select data in specified window
		E_cut, I_cut = self._windowcut(self.E,self.i,E_corr, window)

		# define standard error for fitting
		if apply_weight_distribution:
			sigma = self._apply_weight_distribution(E_cut, I_cut, E_corr, w_ac, W)
		else:
			sigma = np.full((len(E_cut)), 1)

		# Obtain the logarithm for the current densities
		i_corr_guess = math.log10(np.abs(i_corr_guess))
		i_L_guess = math.log10(np.abs(i_L_guess))

		# fitting
		xdata = np.zeros((len(E_cut),2))
		xdata[:,0] = E_cut
		xdata[:,1] = E_corr
		if fix_i_L:
			popt, pcov = curve_fit(forward, xdata, I_cut, sigma=sigma, p0=[i_corr_guess, 0.0600, 0.100, i_L_guess,3 ], bounds=((math.log10(np.abs(self.i).min()),0,0,i_L_guess-0.01,2),(math.log10(np.abs(self.i).max()),1,1,i_L_guess+0.01,4)))
			d_final = forward(xdata, popt[0], popt[1], popt[2], popt[3],popt[4])
			d_final[d_final == 0] = 1E-8  # to get rid of zeros
			RMSE=math.sqrt(np.sum(np.square(I_cut-d_final)).mean())

		else:
			popt, pcov = curve_fit(forward, xdata, I_cut, sigma=sigma, p0=[i_corr_guess, 0.0600, 0.100, i_L_guess,3 ], bounds=((math.log10(np.abs(self.i).min()),0,0,math.log10(np.abs(self.i).min()),2),(math.log10(np.abs(self.i).max()),1,1,math.log10(np.abs(self.i).max())+1,4))) 
			d_final = forward(xdata, popt[0], popt[1], popt[2], popt[3],popt[4])
			d_final[d_final == 0] = 1E-8  # to get rid of zeros
			RMSE=math.sqrt(np.sum(np.square(I_cut-d_final)).mean())
		
		self.fit_results = [d_final, E_cut]
		self.E_corr = E_corr
		self.I_corr = popt[0]
		self.b_a = popt[1]
		self.b_c = popt[2]
		self.i_L = popt[3]
		self.gamma = popt[4]
		self.RMSE = RMSE
		self.w_ac = w_ac
		self.W = W
		self.apply_W = apply_weight_distribution


		if obtain_io:
			io_an, io_cath = self._compute_exchcurrent(popt[0],E_corr, popt[1], E_rev_an, popt[2], E_rev_cath)
			self.io_an = io_an
			self.io_cath = io_cath
			return [d_final, E_cut], E_corr, 10**popt[0], popt[1], -popt[2], 10**popt[3], RMSE, popt[4], io_an, io_cath
		else:
			self.io_an = None
			self.io_cath = None
			return [d_final, E_cut], E_corr, 10**popt[0], popt[1], -popt[2], 10**popt[3], RMSE, popt[4]
	

########################### Sensitivity analysis #############################
##############################################################################

	def sens_analysis(self, window, w_dc, dw_c = 0.01, dw_ac = 0.01, W = np.append(np.arange(50,100,5),[0]), w_c = None, w_ac = None, i_corr_guess = None, i_L_guess = None, fix_i_L = False, output_folder='sensitivity_analysis'):
		
		"""
		Sensitivity analysis of the 'mixed activation-diffusion control fit' to the set parameters, W and w_ac, for the weight distribution. The goal of this function is to find the parameters that produce most stable results vs the cathodic window.
		
		:param window: Lower and upper bounds of the total data to be taken into account in the analysis, relative to the corrosion potential [min value, max value], units [V vs E_corr] 
		:type window: 2-length sequence

		:param w_dc: (An estimation of) the potential versus the corrosion potential (E_corr) at which the diffusion controlled domain starts [V vs E_corr]
		:type w_dc: float

		:param W: array contaning a number M of different values for the weight percentage, W [%], to take into account in the visualization. A W of 0% equals the case where no weight distribution is applied. (Default: [50 55 60 65 70 75 80 85 90 95 0])
		:type W: M-length sequence

		:param w_ac: array contaning a number N of different values for the window activation control, w_ac [V vs OCP], to take into account in the visualization. (Default: 0.01 --> w_dc, increment of 0.01)
		:type w_ac: N-length sequence
		
		:param i_corr_guess: First guess of (the order of magnitude of) the corrosion current density (optional) [A/surface area], it might lead to faster convergence (Default: 10^-2).
		:type i_corr_guess: float
		
		:param i_L_guess: First guess of the limiting current density (optional) [A/surface area], it can to faster convergence and a more accurate fit (Default: 10^1).
		:type i_L_guess: float
		
		:param fix_i_L: If True, the limiting current density is fixed at the value of i_L_guess. (Default:False)
		:type fix_i_L: bool

		:param output_foler: The main output directory in which the the results will be saved (Default:'sensitivity_analysis')
		:type output_folder: string
		
		:returns:
			- fitted_parameters.txt: text file containing the results of the sensitivity analyses, the obtained cathodic Tafel slope and i_L for every evaluated compbination of W, w_ac and window cathodic.
			- stability.txt: The standard deviation of the cathodic Tafel slope data with cathodic window larger than w_dc, for all computed combination of w_ac and W.
			- stability_W.jpeg: Mean standard deviation of the cathodic Tafel slope data with cathodic window larger than w_dc for varying w_ac, as a function of W.
			- stability_wac: Folder containing figures of the standard deviation of the cathodic Tafel slope data with cathodic window larger than w_dc as a function for w_ac, for different W.
		"""

		# Check input values
		if w_dc > 0:
			raise ValueError("w_dc, the start of the diffusion controlled domain should be negative [V vs E_corr]")

		if dw_c < 0:
			raise ValueError("dw_c, the increment of the sweep of the cathodic window should be positive [V]")

		# obtaining automatic guesses for i_icorr and i_L
		if i_corr_guess == None:
			i_corr_guess = np.mean(self.i)
		
		if i_L_guess == None:
			i_L_guess = np.mean(self.i)

		# check input 
		if (i_corr_guess>self.i.max() or i_corr_guess<self.i.min()):
			raise ValueError("Specified i_corr_guess does not lie within the range of the current density")
		if (i_L_guess>self.i.max() or i_L_guess<self.i.min()):
			raise ValueError("Specified i_L_guess does not lie within the range of the current density")

		# making directory
		try:
			os.makedirs(output_folder + '/variability_wac')
		except:
			print("Output folder exists - files are overwritten")

		# Initializing ranges for parameter search
		self.w_dc = w_dc
		W_, w_ac_, window_cat_ = self._param_ranges(W, w_ac, w_c, window, w_dc, dw_c, dw_ac)

		# Initializing panda data frame
		df = pd.DataFrame(data={'w_ac': [],
	                      		'W': [],
	                      		'window_cat': [],
	                      		'b_c': [],
								'i_L': []})

		# Obtaining cathodic tafel slopes, b_c, and limiting current densities, i_L - parameter search
		timer = 0
		for w_ac in w_ac_:
			print('Multi-parameter-sweep in progress, completed:',round((len(w_ac_)-(len(w_ac_)-timer))/len(w_ac_)*100,2),'%       ', end='\r')
			for W in W_:
				for window_cat in window_cat_:
					if -w_ac>window_cat:
						try:
							if W == 0.0:
								result = self.mixed_pol_fit(window=[window_cat,window[1]], i_corr_guess = i_corr_guess, fix_i_L = fix_i_L, i_L_guess = i_L_guess, apply_weight_distribution = False)
							else:
								result = self.mixed_pol_fit(window=[window_cat,window[1]], i_corr_guess = i_corr_guess, fix_i_L = fix_i_L, i_L_guess = i_L_guess, apply_weight_distribution = True, w_ac = w_ac, W = W)
							df_temp = pd.DataFrame(data={'w_ac': [w_ac],
								                         'W': [W],
								                         'window_cat': [window_cat],
								                         'b_c': [result[4]], 
								                         'i_L': [result[5]]})
							df = df.append(df_temp, ignore_index=True)
						except: 
							print("No fit found for W = ", W, ", w_ac = ", w_ac, "and cathodic window = ", window_cat,". Combination skipped.")
			timer += 1
		print('Multi-parameter-sweep in progress, completed:',100.0,'%       ')

		self.df_sens = df
		self.window_sens = window

		# write df to output file
		df.to_csv(output_folder + '/fitted_parameters.txt', sep = '\t', float_format = '%.4g', index=False)


		# Obtain 'variabilit', quantified by the standard deviation of the data for varying window_cat>w_dc
		df_std = pd.DataFrame(data={'w_ac':[],
			                        'W': [],
			                        'std_dev': []})
		for w_ac in w_ac_:
			for W in W_:
				deviation = df.loc[(df['window_cat']<w_dc) & (df['w_ac'] == w_ac) & (df['W'] == W)]['b_c'].std()
				df_std = df_std.append(pd.DataFrame(data={'w_ac':[w_ac],'W': [W],'std_dev': [deviation]}))

		self.df_sens_std = df_std

		# write df_std to output file
		df_std.to_csv(output_folder + '/variability.txt', sep = '\t', float_format = '%.4g', index=False)

		# Plot standard deviation as a function of W and w_ac
		fig,ax = plt.subplots(figsize=(6, 5))
		mean_ = []
		deviation_ = []
		for W in W_:
			if W == 0.0:
				plt.scatter(W,df_std.loc[df_std['W']== W]['std_dev'].mean(), s = 200, c = '#000000', marker = '*')
			else:
				mean_.append(df_std.loc[df_std['W']== W]['std_dev'].mean())
				deviation_.append(df_std.loc[df_std['W']== W]['std_dev'].std())
		W_ = W_[W_ != 0]
		plt.errorbar(W_, mean_,  yerr=deviation_,fmt='.k',  capsize = 2, linewidth=2, elinewidth=1)
		plt.plot(W_, mean_,'-ko')
		plt.xlabel('W [%]')
		plt.ylabel(r'mean variability, $\overline{\sigma^{W}(w_{ac})}$')
		plt.tight_layout()
		fig.savefig(output_folder +'/variability_W.jpeg', format='jpeg', dpi=1000)

		for W in W_:
			plt.close('all')
			fig,ax = plt.subplots(figsize=(6, 5))
			value_array = []
			for w_ac in w_ac_:
			    value_array.append(df_std.loc[(df_std['w_ac']== w_ac) &(df_std['W']== W)]['std_dev'])
			plt.plot(w_ac_, value_array, '-ko')
			plt.xlabel(r'$w_{ac}$ [V]')
			plt.ylabel(r'variability, $\sigma^{W}(w_{ac})$')
			plt.title('W = ' + str(W))
			plt.tight_layout()
			fig.savefig(output_folder +'/variability_wac/W='+'%.1f' % W+'.jpeg', format='jpeg', dpi=1000)

	def plotting_sens_analysis(self, dw_c = 0.01, dw_ac = 0.01,  W = np.append(np.arange(50,100,5),[0]), w_ac = None, w_c = None, output_folder = 'sensitivity_analysis'):
		
		"""
		It returns 3 plots, visualizing the results of the sensitivity analysis: showing the effect of W & w_ac on on i_L and b_cath as a function of the amount of the cathodic branch taken into account in the fitting (cathodic window).
		The 3 plots are safed in different folders in the output_folder
		
		:param W: array contaning a number M of different values for the weight percentage, W [%], to take into account in the visualization. A W of 0% equals the case where no weight distribution is applied. (Default: [50 55 60 65 70 75 80 85 90 95 0])
		:type W: M-length sequence

		:param w_ac: array contaning a number N of different values for the window activation control, w_ac [V vs OCP], to take into account in the visualization. (Default: 0.01 --> w_dc, increment of 0.01)
		:type w_ac: N-length sequence
		
		:param output_folder: The main output directory in which the 5 directories with figures will be saved (Default:'sensitivity_analysis')
		:type output_folder: string

		:return 3 figures:
			- effect_W: The effect of the W on the cathodic Tafel slope b_cath, as a function of the cathodic window (plotted for different w_ac)
			- effect_W_il: The effect of the W on the limiting current density i_L, as a function of the cathodic window (plotted for different w_ac)
			- effect_window_act_control: The effect of w_ac on b_cath, as a function of the cathodic window (plotted for different W)
		"""
		
		# Check input values
		if isinstance(self.df_sens, pd.DataFrame) == False:
			raise ValueError("No results found of the sensitivity analysis.")

		# making the directories
		try:
			os.makedirs(output_folder+'/effect_W_il')
			os.makedirs(output_folder+'/effect_W')
			os.makedirs(output_folder+'/effect_wac')

		except:
			print('Output folder(s) exist(s) - plots will be overwritten')

		# Initializing ranges for parameter search
		W_, w_ac_, window_cat_ = self._param_ranges(W, w_ac, w_c, self.window_sens, self.w_dc, dw_c, dw_ac)

		# plot 1: effect_W
		for w_ac in w_ac_:
			self._plot_effect_W(self.df_sens,w_ac,W_,'b_c',r'$\beta_{cath}$ [V/dec]',output_folder + '/effect_W/wac=')

		# plot 2: effect_W_il
		for w_ac in w_ac_:
			self._plot_effect_W(self.df_sens,w_ac,W_,'i_L',r'$i_{L}$ [A/m$^2$]',output_folder + '/effect_W_il/wac=')

		# plot 3: effect_window_act_control
		for W in W_:
			if W != 0:
				self._plot_effect_Wac(self.df_sens,W,w_ac_,'b_c',r'$\beta_{cath}$ [V/dec]',output_folder+'/effect_wac/W=')
						
############################## Plotting ######################################
##############################################################################

	
	def plotting(self, output_folder = 'plots_out'):

		"""
		Visualization of the results obtained from the fitting. It plots 4 figures:
		data: plot of the raw data as given as input (no IR correction or correction for the surface are applied)
		fit_linear: plot of the observed data and the fit
		fit_semilogarithmic: plot of the observed data and the fit in semilogarithmic scale
		results_overview: plot of the observed data and the obtained Tafel slope(s) (and limiting current density)

		:param output_folder: Directory in which the figures are safed. (Default: 'plots_out')
		:type output_folder: string

		"""

		# Check input value
		if self.b_a == None:
			raise ValueError("No fitting results found to plot.")

		# Make output directory
		try:
			os.makedirs(output_folder)
		except:
			print('Output folder exists - plots will be overwritten')

		# plot 1: data
		plt.close('all')
		fig = plt.figure(1)
		plt.plot(self.E_obs,self.I)
		plt.ylabel('I [A]')
		plt.xlabel('E [V vs Ref]')
		plt.title('Data given as input')
		fig.savefig(output_folder + '/data.jpeg', format='jpeg', dpi=1000)

		# plot 2: fit_linear
		fig = plt.figure(2)
		plt.plot(self.E,self.i, label = 'observed')
		plt.plot(self.fit_results[1],self.fit_results[0], label = 'fit')
		plt.ylabel(r'i [A/m$^2$]')
		plt.xlabel('E [V vs Ref]')
		plt.legend(loc='upper left')
		plt.title('fit')
		fig.savefig(output_folder + '/fit_linear.jpeg', format='jpeg', dpi=1000)

		# plot 3: fit_semilogarithmic
		fig = plt.figure(3)
		plt.plot(self.E,(np.abs(self.i)), label = 'observed')
		plt.plot(self.fit_results[1],np.abs(self.fit_results[0]), label = 'fit')
		plt.ylabel(r'|i| [A/m$^2$]')
		plt.xlabel('E [V vs Ref]')
		plt.yscale('log')
		plt.legend(loc='lower right')
		plt.title('fit')
		fig.savefig(output_folder + '/fit_semilogarithmic.jpeg', format='jpeg', dpi=1000)

		# plot 4: results_overview
		fig = plt.figure(4)
		plt.plot(self.E,np.log10(np.abs(self.i)), '-k', label = 'polarization curve')
		plt.plot(self.E,((self.I_corr)-((self.E_corr)-(self.E))/self.b_a), '-b', label = 'Tafel slope')
		if self.b_c!=None:
			plt.plot(self.E,((self.I_corr)-((self.E_corr)-(self.E))/-self.b_c), '-b') #, label = 'cathodic Tafel slope')
		if self.i_L!= None:
			plt.plot(self.E,([self.i_L]*len(self.E)), '-r', label = 'limiting current density')
		plt.ylabel(r'log10|i| [A/m$^2$]')
		plt.xlabel('E [V vs Ref]')
		plt.legend(loc='lower left')
		plt.title('Overview results')
		fig.savefig(output_folder + '/results_overview.jpeg', format='jpeg', dpi=1000)
		plt.close('all')
		
############################## Writing ######################################
##############################################################################

	def save_to_txt(self, filename = './fitting_results'):

		"""
		Saving of the results obtained from the fitting in a text file. It lists the fitting technique used, the fitted parameters, and 
		the corresponding fitted curve. 

		:param filename: path and name of the output file (Default: './fitting_results')
		:type output_folder: string

		"""
		
		# Check input data
		if self.b_a == None:
			raise ValueError("No fitting results found.")

		with open(filename + '.txt', 'w') as f:
			if self.b_c == None:
				f.write('Results, obtained by the linear fit. \n\n Fitted parameters: ' + '\n')
			elif self.i_L == None:
				f.write('Results, obtained by the activation control fit \n\n Fitted parameters: ' + '\n')
			else:
				f.write('Results, obtained by the mixed diffusion-activation control fit \n\n')
				if self.apply_W:
					f.write('Weight distribution applied: \n')
					f.write('window activation control (w_ac) [V] = '+ '%.3f' % self.w_ac + '\n')
					f.write('Weight percentage (W) [%] = '+ '%.1f' % self.W + '\n\n')

				f.write('Fitted parameters: ' + '\n')
				
			f.write('Corrosion potential (E_corr) [V vs Ref] = ' + '%.4f' % self.E_corr + '\n') 
			f.write('Corrosion current (density) (I_corr)  = ' + '%.4f' % 10**self.I_corr + '\n')
			if self.b_c == None:
				f.write('Tafel slope (b) [V]  = ' + '%.4f' % self.b_a + '\n')
			else:
				f.write('Anodic Tafel slope (b_a) [V] = ' + '%.4f' % self.b_a + '\n')
				f.write('Cathodic Tafel slope (b_c) [V] = ' + '%.4f' % -self.b_c + '\n')
			if self.i_L != None:
				f.write('Limiting current (density) (i_L) = ' + '%.4f' % 10**self.i_L + '\n')
			f.write('RMSE = ' + '%.4f' % self.RMSE + '\n')
			if self.b_c == None:
				if self.io_an != None:
					f.write('Exchange current (density) (io) = ' + '%.4f' % self.io_an + '\n')
			else:
				if self.io_an != None:
					f.write('Anodic exchange current (density) (io_an) = ' + '%.4f' % self.io_an + '\n')
				if self.io_cath != None:
					f.write('Cathodic exchange current (density) (io_cath) = ' + '%.4f' % self.io_cath + '\n')
			f.write('\n')
			f.write('Fitted curve: \n')
			f.write('E_fit [V vs Ref] \t I_fit [Unit input current/ input surface area] \n')
			np.savetxt(f, np.c_[self.fit_results[1], self.fit_results[0]],newline='\n')



	############################## Data functions ################################
	##############################################################################

	# function to obtain the weight distribution for the 'mixed activation-diffusion control fit'
	def _apply_weight_distribution(self, E, I, E_corr, w_ac, W):
		E_sigma, I_sigma = self._windowcut(E, I, E_corr,[-w_ac,w_ac])
		sigma = np.full(len(E), 1.0)
		for j in range(len(E)):
			if (-w_ac <= (E[j]-E_corr)) and ((E[j]-E_corr)<= w_ac):
				sigma[j] = (100-W)/W * len(E_sigma)/(len(E)-len(E_sigma))

		return sigma

	# function to correct the data for the voltage drop (IR correction)
	def _ir_correction(self,R):
		return self.E-self.I*R

	# function to obtain the corrosion potential E_corr
	def _find_Ecorr(self):
		min_index = np.argmin(np.abs(self.i)) 
		return self.E[min_index]

	# function to select the data within a certain window
	def _windowcut(self, E, i, E_corr, window):
		E_cut = []
		I_cut = []
		for j in range(len(i)):
			if (window[0] <= (E[j]-E_corr)) and ((E[j]-E_corr)<= window[1]):
				E_cut.append(E[j])
				I_cut.append(i[j])
		return E_cut, I_cut

	# function to compute the exchange current densities
	def _compute_exchcurrent(self,I_corr, E_corr, b_a, E_rev_an, b_c = None, E_rev_cath = None):
	    if b_c == None:
	    	io = I_corr - 1/b_a*(E_rev_an-0) # 0 corresponding to the computation of I_corr
	    	return io
	    else:
	    	io_an = I_corr + 1/b_a*(E_rev_an-E_corr)
	    	io_cath = I_corr - 1/b_c*(E_rev_cath-E_corr)
	    return io_an, io_cath

	# function to obtain the change in respect to the mean of the last 10 points
	def _get_fluc(self,y):
		y_fluc = np.zeros(len(y))
		for j in range(1,len(y)-1):
			if j<9:
				y_fluc[j+1] = (y[j+1]-np.mean(y[:j]))
			else:
				y_fluc[j+1] = (y[j+1]-np.mean(y[j-9:j]))
		return y_fluc
				
	# function to truncate a colormap
	def _truncate_colormap(self, cmap, minval=0.0, maxval=1.0, n=100):
		new_cmap = colours.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),cmap(np.linspace(minval, maxval, n)))
		return new_cmap

	# function for plotting, dependency of W
	def _plot_effect_W(self,df,w_ac,W_,param,ylabel,output_folder):
		plt.close('all')
		fig,ax = plt.subplots(figsize=(10, 5))
		for W in W_:				
			df_select = df.loc[df['w_ac'] == w_ac]
			df_select = df_select.loc[df_select['W']== W]
			x = df_select['window_cat'].to_numpy()
			y = df_select[param].to_numpy()

			if W == 0:
				plt.plot(x,y,'-k', label = 'not weighted')
			else:
				plt.plot(x,y, label = 'W = '+ str(W) +' %')

			if param == 'i_L':
				plt.yscale('log')

		plt.ylabel(ylabel)
		plt.xlabel(r'cathodic window $w_{c}$ [V from $E_{corr}$]')
		plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
		plt.title(r'window activation control $w_{ac}$ = ' + '%.2f' % w_ac + ' V')
		plt.tight_layout()
		fig.savefig(output_folder +'%.2f' % w_ac+'.jpeg', format='jpeg', dpi=1000)

	# function for plotting, dependency on w_ac
	def _plot_effect_Wac(self,df,W,w_ac_,param,ylabel,output_folder):
		# defining colorscale 
		colors = plt.cm.Reds(np.linspace(0.2, 1.0, len(w_ac_)))
		plt.close('all')
		fig,ax = plt.subplots(figsize=(10, 6))
		i = 0
		for w_ac in w_ac_:
			df_select = df.loc[df['w_ac'] == w_ac]	
			df_select = df_select.loc[df_select['W'] == W]
			x = df_select['window_cat'].to_numpy()
			y = df_select[param].to_numpy()
			
			plt.plot(x,y, color=colors[i], label = str(w_ac))
			i+=1
		df_select_notweighted = df.loc[(df['W']==0) & (df['w_ac']==0.01)]
		plt.plot(df_select_notweighted['window_cat'].to_numpy(),df_select_notweighted[param].to_numpy(), '--k')
		cmap=self._truncate_colormap(plt.get_cmap('Reds'), 0.2, 1.0)
		norm = mpl.colors.Normalize(vmin=w_ac_.min(),vmax=w_ac_.max())
		plt.ylabel(ylabel)
		plt.xlabel(r'cathodic window $w_{c}$ [V from $E_{corr}$]')
		plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),label=r'w$_{ac}$ [V]')
		plt.title('Weight (W) = ' + str(W)+'%')
		fig.savefig(output_folder+str(W)+'.jpeg', format='jpeg', dpi=1000)


	# function to define and check ranges for parameter search
	def _param_ranges(self, W, w_ac, w_c, window, w_dc, dw_c, dw_ac):
		
		# Initializing ranges parameter search
		W_ = np.array(W)

		w_ac_ = np.array(w_ac)
		try:
			if w_ac == None:
				w_ac_ = np.arange(dw_ac,round(abs(0.5*w_dc)+dw_ac,2),dw_ac)
		except:
			pass
		
		w_c_ = np.array(w_c)
		try:
			if w_c == None:
				w_c_ = np.arange(w_dc, round(window[0]-dw_c,2), -dw_c)
		except:
			pass

		# Check input values and ranges
		if np.any(( W_< 0)|(W_ > 100)):
			raise ValueError("Values for W should be in the range of  0-100 %")
		if np.any(w_ac_ < 0):
			raise ValueError("Values for w_ac should be larger than 0")

		return W_, w_ac_, w_c_
	
	


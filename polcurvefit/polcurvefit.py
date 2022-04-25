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
import warnings

class PolCurveFit:
	
	"""
	PolCurveFit is a python class that can be used to analyze measured polarization curves and obtain parameters such as the corrosion potential, 
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

		# Initializing fitting parameters
		self.fit_results = None
		self.E_corr = None
		self.I_corr = None
		self.b_a = None
		self.b_c = None
		self.i_L = None
		self.gamma = None
		self.RMSE = None
		self.io_an = None
		self.io_cath = None
		self.w_ac = None
		self.W = None
		self.apply_W = None

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
		if (i_corr_guess>self.i.max() or i_corr_guess<self.i.min()) and i_corr_guess!= 10**-2.103:
			warnings.warn("Specified i_corr_guess does not lie within the range of the current density")

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
		if (i_corr_guess>self.i.max() or i_corr_guess<self.i.min()) and i_corr_guess!= 10**-2.103:
			warnings.warn("Specified i_corr_guess does not lie within the range of the current density")
		if (i_L_guess>self.i.max() or i_L_guess<self.i.min()) and i_L_guess!= 10**1.103:
			warnings.warn("Specified i_L_guess does not lie within the range of the current density")

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

	def sens_analysis(self, window, start_diffusion_control, importance = [50,70,80,90,0], i_corr_guess = 10**-2, i_L_guess = 10**1, fix_i_L = False, output_folder='sensitivity_analysis'):
		
		"""
		Sensitivity analysis of the 'mixed activation-diffusion control fit' to the set parameters for the weight distribution.
		It returns 5 plots, showing the effect of the importance & w_acon on i_L and b_cath as a function of the amount of the cathodic branch taken into account in the fitting (cathodic window).
		The 5 plots are safed in different folders in the output_folder

		
		:param window: Lower and upper bounds of the total data to be taken into account in the analysis, relative to the corrosion potential [min value, max value], units [V vs E_corr] 
		:type window: 2-length sequence

		:param start_diffusion_control: (An estimation of)the potential versus the corrosion potential (E_corr) at which the diffusion controlled domain starts [V vs E_corr]
		:type start_diffusion_control: float

		:param importance: array contaning a number M of different values for importance [%] to take into account in the analysis. An importance of 0% equals the case where no weight distribution is applied. (Default: [50,70,80,90,0])
		:type importance: M-length sequence
		
		:param i_corr_guess: First guess of (the order of magnitude of) the corrosion current density (optional) [A/surface area], it might lead to faster convergence (Default: 10^-2).
		:type i_corr_guess: float
		
		:param i_L_guess: First guess of the limiting current density (optional) [A/surface area], it can to faster convergence and a more accurate fit (Default: 10^1).
		:type i_L_guess: float
		
		:param fix_i_L: If True, the limiting current density is fixed at the value of i_L_guess. (Default:False)
		:type fix_i_L: bool

		:param output_foler: The main output directory in which the 5 directories with figures will be saved (Default:'sensitivity_analysis')
		:type output_folder: string

		:return 5 figures:
			- effect_importance: The effect of the importance on the cathodic Tafel slope b_cath, as a function of the cathodic window (plotted for different w_ac)
			- effect_importance_fluctuation: The effect of the importance on the change in b_cath in respect to the mean of the previous 100 mV of smaller absolute cathodic windows (plotted for different w_ac)
			- effect_importance_il: The effect of the importance on the limiting current density i_L, as a function of the cathodic window (plotted for different w_ac)
			- effect_window_act_control: The effect of w_ac on b_cath, as a function of the cathodic window (plotted for different importance)
			- effect_window_act_control_fluctuation: The effect of w_ac on the change in b_cath in respect to the mean of the previous 100 mV of smaller absolute cathodic windows (plotted for different importance)
		"""

		# Initializing parameter search
		window_cat_ = np.arange(-0.020, window[0]-0.010, -0.010)
		w_ac_ = np.arange(0.010,abs(start_diffusion_control)+0.010,0.010)
		importance_ = importance

		# Initializing panda data frame
		df = pd.DataFrame(data={'b_c_best': [],
								'i_L_best': [],
	                       		'window_cat':[], 
	                      		'w_ac': [],
	                      		'importance':[]})

		# Obtaining fitting results - parameter search
		for window_cat in window_cat_:
			for w_ac in w_ac_:
				for importance in importance_:
					if -w_ac>window_cat:
						if importance == 0:
							temp_result = self.mixed_pol_fit(window=[window_cat,window[1]], i_corr_guess = i_corr_guess, fix_i_L = fix_i_L, i_L_guess = i_L_guess, apply_weight_distribution = False)
						else:
							temp_result = self.mixed_pol_fit(window=[window_cat,window[1]], i_corr_guess = i_corr_guess, fix_i_L = fix_i_L, i_L_guess = i_L_guess, apply_weight_distribution = True, w_ac = w_ac, W = importance)
						df_temp = pd.DataFrame(data={'b_c_best': [temp_result[4]], 
							                         'i_L_best':[temp_result[5]],
							                         'window_cat':[window_cat], 
							                         'w_ac': [w_ac],
							                         'importance':[importance]})
						df = df.append(df_temp, ignore_index=True)

		# plotting

		# making the directories
		try:
			os.makedirs(output_folder+'/effect_importance_il')
			os.makedirs(output_folder+'/effect_importance')
			os.makedirs(output_folder+'/effect_importance_fluctuation')
			os.makedirs(output_folder+'/effect_window_act_control')
			os.makedirs(output_folder+'/effect_window_act_control_fluctuation')

		except:
			print('Warning: Output folder(s) exist(s) - plots will be overwritten')
		
		# plot 1: effect_importance
		for w_ac in w_ac_:
			self._plot_effect_W(df,w_ac,importance_,'b_c_best',r'$\beta_{cath}$ [V/dec]',output_folder + '/effect_importance/wac=',fluctuation=False)

		# plot 2: effect_importance_fluctuation
		for w_ac in w_ac_:
			self._plot_effect_W(df,w_ac,importance_,'b_c_best',r'relative change of $\beta_{cath}$ [V/dec]',output_folder + '/effect_importance_fluctuation/wac=',fluctuation=True)

		# plot 3: effect_importance_il
		for w_ac in w_ac_:
			self._plot_effect_W(df,w_ac,importance_,'i_L_best',r'$i_{L}$ [A/m$^2$]',output_folder + '/effect_importance_il/wac=',fluctuation=False)

		# plot 4: effect_window_act_control
		for importance in importance_:
			if importance != 0:
				self.__plot_effect_Wac(df,importance,w_ac_,'b_c_best',r'$\beta_{cath}$ [V/dec]',output_folder+'/effect_window_act_control/W=',fluctuation=False)
			
		# plot 5: effect_window_act_control_fluctuation
		for importance in importance_:
			if importance != 0:
				self.__plot_effect_Wac(df,importance,w_ac_,'b_c_best',r'relative change of $\beta_{cath}$ [V/dec]',output_folder+'/effect_window_act_control_fluctuation/W=',fluctuation=True)
			
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
		the corresponding fitted curve

		:param filename: path and name of the output file (Default: './fitting_results')
		:type output_folder: string

		"""

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
			f.write('Corrosion current (density) (I_corr)  = ' + '%.4f' % self.I_corr + '\n')
			if self.b_c == None:
				f.write('Tafel slope (b) [V]  = ' + '%.4f' % self.b_a + '\n')
			else:
				f.write('Anodic Tafel slope (b_a) [V] = ' + '%.4f' % self.b_a + '\n')
				f.write('Cathodic Tafel slope (b_c) [V] = ' + '%.4f' % -self.b_c + '\n')
			if self.i_L != None:
				f.write('Limiting current (density) (i_L) = ' + '%.4f' % self.i_L + '\n')
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
	def _plot_effect_W(self,df,w_ac,importance_,param,ylabel,output_folder,fluctuation=False):
		plt.close('all')
		fig,ax = plt.subplots(figsize=(10, 5))
		for importance in importance_:				
			df_select = df.loc[df['w_ac'] == w_ac]
			df_select = df_select.loc[df_select['importance']== importance]
			x = df_select['window_cat'].to_numpy()
			y = df_select[param].to_numpy()

			if fluctuation:	
				y = self._get_fluc(y)

			if importance == 0:
				plt.plot(x,y,'-k', label = 'not weighted')
			else:
				plt.plot(x,y, label = 'W = '+ str(importance) +' %')

			if param == 'i_L_best':
				plt.yscale('log')
		plt.ylabel(ylabel)
		plt.xlabel('cathodic window [V from OCP]')
		plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
		plt.title('window activation control = ' + '%.2f' % w_ac + ' mV')
		plt.tight_layout()
		fig.savefig(output_folder +'%.2f' % w_ac+'.jpeg', format='jpeg', dpi=1000)

	# function for plotting, dependency on w_ac
	def _plot_effect_Wac(self,df,importance,w_ac_,param,ylabel,output_folder,fluctuation=False):
		# defining colorscale 
		colors = plt.cm.Reds(np.linspace(0.2, 1.0, len(w_ac_)))
		plt.close('all')
		fig,ax = plt.subplots(figsize=(10, 6))
		i = 0
		for w_ac in w_ac_:
			df_select = df.loc[df['w_ac'] == w_ac]	
			df_select = df_select.loc[df_select['importance'] == importance]
			x = df_select['window_cat'].to_numpy()
			y = df_select[param].to_numpy()
			
			if fluctuation:	
				y = self._get_fluc(y)

			plt.plot(x,y, color=colors[i], label = str(w_ac))
			i+=1

		df_select_notweighted = df.loc[(df['importance']==0) & (df['w_ac']==0.01)]
		if fluctuation:
			plt.plot(df_select_notweighted['window_cat'].to_numpy(),self._get_fluc(df_select_notweighted[param].to_numpy()), '--k')
		else:
			plt.plot(df_select_notweighted['window_cat'].to_numpy(),df_select_notweighted[param].to_numpy(), '--k')
		cmap=self._truncate_colormap(plt.get_cmap('Reds'), 0.2, 1.0)
		norm = mpl.colors.Normalize(vmin=w_ac_.min(),vmax=w_ac_.max())
		plt.ylabel(ylabel)
		plt.xlabel('cathodic window [V from OCP]')
		plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),label=r'w$_{ac}$ [V]')
		plt.title('Weight (W) = ' + str(importance)+'%')
		fig.savefig(output_folder+str(importance)+'.jpeg', format='jpeg', dpi=1000)


	

	
	


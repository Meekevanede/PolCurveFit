import numpy as np
import warnings

#suppress warnings
warnings.filterwarnings('ignore')

# with concentration limited current
def forward( xdata, I_corr, b_a, b_c, i_L,gamma):
	E, E_corr = xdata[:,0], xdata[1,1]
	solution = 10**I_corr * np.exp(2.3*(E-E_corr)/b_a) - (((10**I_corr*np.exp(2.3*(E_corr-E)/b_c))**gamma)/(1+((10**I_corr*np.exp(2.3*(E_corr-E)/b_c))/10**i_L)**gamma))**(1/gamma)
	return solution

# without concentration limited current
def forward_normal( xdata, I_corr, b_a, b_c):
	E, E_corr = xdata[:,0], xdata[1,1]
	solution = 10**I_corr * np.exp(2.3*(E-E_corr)/b_a) - 10**I_corr * np.exp(2.3*(E_corr-E)/b_c)
	return solution

# linear fit
def linear_fit(xdata, b_c, Constant):
	E = np.array(xdata)
	solution =  b_c*E + Constant 
	return solution
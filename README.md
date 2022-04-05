# PolCurveFit
A python library to analyse polarization curves, by fitting theoretical curves to input data. Parameters such as the corrosion potential, corrosion rate, Tafel slopes and exchange current densities can be obtained, with three included techniques:
Tafel extrapolation: a linear fit to a defined Tafel region
Activation control fit: fitting of a theoretical curve describing the anodic and cathodeic activation controlled currents around OCP.
Mixed activation-diffusion control fit: fitting of a theoretical curve describing an anodic domain with solely activation controlled currents and a cathodic domain with (mixed) activation and diffusion controlled currents

### Installation
At this moment (testing phase): copy the package (folder polcurvefit), run the example code below in the same folder as where the package (folder polcurvefit) is located. 

When uploaded on Pypi, it can be installed as follows:

```
pip install pol-curve-fit
```

### Get started
Example of how to apply the code

```Python
from polcurvefit import PolCurveFit
import numpy as np

# upload an example polarization curve
input_data = np.loadtxt('polcurvefit/data/example_data.txt',comments = '#', delimiter = '\t', usecols=(0, 1))
E=input_data[:,0]
I=input_data[:,1]

# Instantiate a polarization curve object
Polcurve = PolCurveFit(E,I,sample_surface=2.0106E-04)

# Apply a fitting technique: 'the activation control fit':
fitted_curve, E_corr, I_corr, b_a, b_c, RMSE = Polcurve.active_pol_fit(window=[-0.1,0.04], i_corr_guess = 10**-2)

# Visualise the obtained fit
Polcurve.plotting(output_folder='Visualization_activation_fit')

# Apply a fitting technique: 'the mixed activation-diffusion control fit':
fitted_curve, E_corr, I_corr, b_a, b_c, i_L, RMSE, gamma = Polcurve.mixed_pol_fit(window=[-0.4,0.04], i_corr_guess = 10**-2, apply_weight_distribution = True, w_ac = 0.04, W = 85)

# Visualise the obtained fit
Polcurve.plotting(output_folder='Visualization_mixed_fit')

```
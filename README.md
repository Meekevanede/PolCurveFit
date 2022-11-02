# PolCurveFit
A python library to analyse polarization curves, by fitting theoretical curves to input data. Parameters such as the corrosion potential, corrosion rate, Tafel slopes and exchange current densities can be obtained, with three included techniques:
Tafel extrapolation: a linear fit to a defined Tafel region
Activation control fit: fitting of a theoretical curve describing the anodic and cathodeic activation controlled currents around OCP.
Mixed activation-diffusion control fit: fitting of a theoretical curve describing an anodic domain with solely activation controlled currents and a cathodic domain with (mixed) activation and diffusion controlled currents

### Installation

```
pip install PolCurveFit
```

### Documentation
Find the documentation on: https://polcurvefit.readthedocs.io/

### Example
Example of how to apply the code

```Python
import numpy as np

# An artificial polarization curve
E = np.arange(-1.0, 0.1, 0.01)
I = 0.002 * np.exp(2.3*(E+0.2)/0.08) - 0.002*np.exp(2.3*(-0.2-E)/0.18)/(1+((0.002*np.exp(2.3*(-0.2-E)/0.18))/0.3))

# Instantiate a polarization curve object
from polcurvefit import polcurvefit
Polcurve = polcurvefit(E,I, R= 0, sample_surface=1E-3)

# Apply a fitting technique: 'the activation control fit':
results = Polcurve.active_pol_fit(window=[-0.3,0.1])

# Save the results and visualise the obtained fit
Polcurve.save_to_txt(filename = './output_act')
Polcurve.plotting(output_folder='figures/output_act')

# Apply a fitting technique: 'the mixed activation-diffusion control fit' with a specific weight distribution:
results = Polcurve.mixed_pol_fit(window=[-0.7,0.1], apply_weight_distribution = True, w_ac = 0.07, W = 80)

# Save the results and visualise the obtained fit
Polcurve.save_to_txt(filename = './output_mixed')
Polcurve.plotting(output_folder='figures/output_mixed')

```

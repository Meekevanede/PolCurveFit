# PolCurveFit
A python library to analyse polarization curves, by fitting theoretical curves to input data. Parameters such as the corrosion potential, corrosion rate, Tafel slopes and exchange current densities can be obtained, with three included techniques:
Tafel extrapolation: a linear fit to a defined Tafel region
Activation control fit: fitting of a theoretical curve describing the anodic and cathodeic activation controlled currents around OCP.
Mixed activation-diffusion control fit: fitting of a theoretical curve describing an anodic domain with solely activation controlled currents and a cathodic domain with (mixed) activation and diffusion controlled currents

### Installation
```
pip install pol-curve-fit
```

### Get started
Example of how to apply the code

```Python
from polcurvefit import PolCurveFit

# upload an example polarization curve
E,I = example

# Instantiate a polarization curve object
Polcurve = PolCurveFit(E,I)

# Apply a fitting technique
fitted_curve, E_corr, I_corr, b_a, b_c, RMSE = Polcurve.active_pol_fit(window=[-0.4,0.04], i_corr_guess = 10*-2)

# Visualise the obtained fit
Polcurve.plotting(output_folder='Visualization')
```
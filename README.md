# PolCurveFit
A python library to analyse polarization curves, by fitting theoretical curves to input data. Parameters such as the corrosion potential, corrosion rate, Tafel slopes and exchange current densities can be obtained, with three included techniques:
Tafel extrapolation: a linear fit to a defined Tafel region
Activation control fit: fitting of a theoretical curve describing the anodic and cathodeic activation controlled currents around OCP.
Mixed activation-diffusion control fit: fitting of a theoretical curve describing an anodic domain with solely activation controlled currents and a cathodic domain with (mixed) activation and diffusion controlled currents

### Installation
At this moment (testing phase): copy the package (folder polcurvefit), run the example code below in the same folder as where the package (folder polcurvefit) is located. 

When uploaded on Pypi, it can be installed as follows:

```
pip install PolCurveFit
```

### Documentation
Find the documentation on: https://polcurvefit.readthedocs.io/

### Example
Example of how to apply the code

```Python
import numpy as np

# upload an example polarization curve, using the module DataImport
from polcurvefit import DataImport as di
inputfile = './polcurvefit/data/example_txt.txt'
E,I = di.load_txt(inputfile, lines_header = 1)

# Instantiate a polarization curve object
from polcurvefit import PolCurveFit
Polcurve = PolCurveFit(E,I, R= 0, sample_surface=2.0106E-04)

# Apply a fitting technique: 'the activation control fit':
results = Polcurve.active_pol_fit(window=[-0.05,0])

# Save the results and visualise the obtained fit
Polcurve.save_to_txt(filename = './results/output_act)
Polcurve.plotting(output_folder='figures/output_act')

# Apply a fitting technique: 'the mixed activation-diffusion control fit' with a specific weight distribution:
results = Polcurve.mixed_pol_fit(window=[-0.4,0.1], apply_weight_distribution = True, w_ac = 0.07, W = 80)

# Save the results and visualise the obtained fit
Polcurve.save_to_txt(filename = './results/output_mixed)
Polcurve.plotting(output_folder='figures/output_mixed')

```
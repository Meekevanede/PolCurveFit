.. PolCurveFit documentation master file, created by
   sphinx-quickstart on Tue Apr  5 17:39:40 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to polcurvefit's documentation!
=======================================

The python library polcurvefit allows the user to analyse measured polarization curves and obtain parameters such as the Tafel slopes, exchange current densities and corrosion rates. It implements three techniques: Tafel extrapolation (linear fit), 'activation control fit' and 'mixed activation-diffusion control fit'[van Ede & Angst, 2022]. Here you will find information about the working of the different techniques and how to use the library.

*van Ede, Meeke, and Ueli Angst. "Analysis of polarization curves under mixed activation-diffusion control–an algorithm to minimize human factors." CORROSION (2022): 4171.* 'https://doi.org/10.5006/4171 <https://doi.org/10.5006/4171>'_ *

.. toctree::
   :maxdepth: 2
   :caption: Python Class PolCurveFit

   source/api/introduction
   source/api/polcurvefit
   source/api/dataimport
   
   

.. toctree::
   :maxdepth: 2
   :caption: Getting started
	
   source/examples/getting_started

.. toctree::
   :maxdepth: 2
   :caption: Examples

   source/examples/example1
   source/examples/example2
   source/examples/example3



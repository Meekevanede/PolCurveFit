==============
Methodology
==============

PolCurveFit implements three different fitting techniques: Tafel extrapolation (linear fit), activation control fit & mixed activation-diffusion control fit. These techniques fit a theoretical curve to the polarization curve data. Below a short overview overview of the used theoretical curves for each technique.

Initialisation and correction of the polarization curve
=======================================================

When the polarization curve object is initialised, the user has the possibility to correct the data. Firstly, by giving the electrolyte resistance R as an input, the potentials (E) are corrected for the IR drop:

.. math::

   E = E - I*R

Where I is the current. Secondly the user can specify the surface area of the sample, to be able to convert from currents to current densities. If the input data consists of current densities already, simply use the default value of 1.

Tafel extrapolation (linear fit)
================================
The Tafel extrapolation technique fits a linear line to the polarization curve data within a specified window. The technique is suitable for data that is solely under activation control on the  anodic or on the cathodic branch. This means that only one branch (the cathodic or anodic) can be fitted at a time. If The corrosion potential is given as an input, the fitted Tafel line is extrapolated to also obtain the corrosion current (density).

.. image:: linear.jpeg
   :width: 400
   :alt: Alternative text

Activation control fit
======================
This technique is suitable to fit currents around the corrosion potential which are under activation control. It fits the following equation to the polarization curve data:

.. math::
   
   	i = i_{corr}\left[ \exp\left(\frac{E-E_{corr}}{\beta_{an}}\right)  - \exp\left(\frac{E_{corr}-E}{\beta_{cath}}\right) \right]

.. image:: activation.jpeg
   :width: 400
   :alt: Alternative text

Mixed activation-diffusion control fit
======================================
This technique is suitable for polarization curve data of which the cathodic branch is under mixed activation-diffusion control. It fits the following equation to the data:

.. math::

   i = i_{corr} * \exp{\frac{2.303(E-E_{corr})}{\beta_{an}}} -  \left[ \frac{\left( i_{corr}\exp{\frac{2.303(E_{corr}-E)}{\beta_{cath}}}\right)^\gamma}{1+\left(\frac{i_{corr}}{i_{L}}\exp{\frac{2.303(E_{corr}-E)}{\beta_{cath}}}\right)^\gamma}              \right]^{\frac{1}{\gamma}}

.. image:: mixed.jpeg
   :width: 400
   :alt: Alternative text

The user also has the obtain the use a specific weight distribution. Explain problem and solution, and explain the parameters. --> Use of sensitivity study.
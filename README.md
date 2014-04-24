GaussianCharges
===============

It can be shown that a careful modification of the PME method for point charge coulomb interactions in periodic boundary onditions can yield an equivalent computation for charges distributed as gaussians, thus preserving the O(NlogN) scaling of PME.

Based on that observation, I develop a water model with gaussian charge interactions and the computational efficieny of the PME method in standard models like TIP3P and TIP4P. 

The development is done with the OpenMM  toolkit for molecular simulation, which can be downloaded at https://simtk.org/home/openmm

The water model is still going through  preliminary testing. Once validated, a basic implementation will become available on here.

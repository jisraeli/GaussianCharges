GaussianCharges
===============

It can be shown that a careful modification of the PME method for point charge coulomb interactions in periodic boundary onditions can yield an equivalent computation for charges distributed as gaussians, thus preserving the O(NlogN) scaling of PME.

Based on that observation, I develop a water model with gaussian charge interactions and the computational efficieny of the PME method in standard models like TIP3P and TIP4P. 

The development is done with the OpenMM  toolkit for molecular simulation, which can be downloaded at https://simtk.org/home/openmm

GaussianPME.py simulates water with gaussian charges using OpenMM's Reference and CPU platforms. Note: on a linux machine it may run as much as 100 times slower than an equivalent water model using PME on the CPU platform.

FastGaussianPME.py is an optimized version of GaussianPME.py and it can be implemented on any platform. For large water systems (>1000 particles), it can run as much as 100 times faster than GaussianPME.py

GaussianCharges
===============

Most molecular dynamics simulations use the Particle-Mesh Ewald (PME) method to evaluate coulomb interactions in periodic boundary conditions. The O(NlogN) scaling of the PME algorithm has made it possible to simulate large systems with standard water models such as TIP3P. In this project, I develop a forcefield that models charges as gaussian charges instead of point charges: the goal is a more accurate water model than TIP3P with similar computational efficiency. 

It is possible to modify the PME algorithm so as to compute the gaussian electric potential and preserve the O(NlogN) scaling. Based on that observation, I've developed code using the OpenMM toolkit for molecular simulation which simulates water with gaussian charges. FastGaussianPME.py is the most optimized implementation currently and runs at most x1.2 slower than TIP3P with PME on Linux machines using the OpenCL platform.

temp.xml is a template forcefield xml file that can replace water models in OpenMM and will run at most x1.2 slower than TIP3P with PME. I'm currently optimizing the parameters in temp.xml using the ForceBalance software for systmatic forcefield optimization, you can follow my porgress here: https://github.com/jisraeli/forcebalance

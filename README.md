GaussianCharges
===============

Most molecular dynamics simulations use the Particle-Mesh Ewald (PME) method to evaluate coulomb interactions in periodic boundary conditions. The O(NlogN) scaling of the PME algorithm has made it possible to simulate large systems with standard water models such as TIP3P. In this project, I develop a TIP3G forcefield that models point charges as gaussian charges: the goal is a more accurate water model than TIP3P with similar computational efficiency. 

It is possible to modify the PME algorithm so as to compute the gaussian electric potential and preserve the O(NlogN) scaling. Based on that observation, I've developed code using the OpenMM toolkit for molecular simulation which simulates water with gaussian charges as fast as TIP3P with the PME method (at most 1.2x slower as far as I've checked).

I'm currently in the process of optimizing the forcefield using the ForceBalance software for systematic forcefield optimization. In the results folder you will find reports that include comparisons between my predictions and experimental measurements of 6 standard water properties. The 3Param report is for a version of TIP3G where 3 parameters are being optimized and the 6Param report is for a version where 6 parameters are being optimized. Both reports will be updated as the forcefields improve - the most recent comparisons can be seen at the end of each report.

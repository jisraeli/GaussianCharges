import pandas as pd
import numpy as np
from numpy import average as avg
from math import*
from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import*
from simtk import unit

epsilon = 8.854187817620E-12*farad/meter
COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)
MU_GAS = 1.854989
ALPHA_GAS = 1.6633E-40
kB = BOLTZMANN_CONSTANT_kB 

InFile = './examples/WaterBox.pdb'
pdb = app.PDBFile(InFile)
pdb.topology.setUnitCellDimensions((2.48,2.48,2.48))
forcefield = app.ForceField('tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

N_PARTICLES = system.getNumParticles()

data = pd.read_csv('Reference.csv', sep='\t')
dataGaussian = pd.read_csv('Gaussian.csv', sep='\t')
cols = list(data.columns.values)
density = data[cols[6]]
densityGaussian = dataGaussian[cols[6]]
volume = N_PARTICLES/density
volumeGaussian = N_PARTICLES/densityGaussian
volume2 = volume**2
volume2Gaussian = volumeGaussian**2
PE = data[cols[2]]
PE_Gaussian = dataGaussian[cols[2]]
temp = data[cols[5]]
tempGaussian= dataGaussian[cols[5]]

print "TIP3P results: "
print "\tAverage Density: ", avg(density)
print "\tAverage PE per molecule: ", avg(PE)/(N_PARTICLES/3.0)
print "\tIsothermal Compressibility: ", (avg(volume2)-avg(volume)**2)/(kB*avg(temp)*avg(volume))

print "GaussianCharges results:"
print "\tAverage Gaussian Density: ", avg(densityGaussian)
print "\tAverage PE per molecule: ", avg(PE_Gaussian)/(N_PARTICLES/3.0)
print "\tIsothermal Compressibility: ", (avg(volume2Gaussian)-avg(volumeGaussian)**2)/(kB*avg(tempGaussian)*avg(volumeGaussian)) 
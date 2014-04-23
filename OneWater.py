from math import*
import numpy as np
from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import*
from simtk import unit
from pylab import*
from sys import stdout
import simulation1
import sys

pdb = app.PDBFile('OneWater.pdb')
pdb.topology.setUnitCellDimensions((2,2,2))
forcefield = app.ForceField('tip3p.xml')

systemTest = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

PME = systemTest.getForce(2)
N_PARTICLES = systemTest.getNumParticles()
print "#of particles: ", N_PARTICLES
for i in range(N_PARTICLES):
    charge ,sigma, epsilon = PME.getParticleParameters(i)
    PME.setParticleParameters(i, charge, sigma, 0.0)

N_EXCEPTIONS = PME.getNumExceptions()
print "#of exceptions: ", N_EXCEPTIONS
print "exception parameters: "
for i in range(N_EXCEPTIONS):
    print PME.getExceptionParameters(i)

integratorTest = mm.VerletIntegrator(0.002)
for i in range(len(systemTest.getForces())):
    force = systemTest.getForce(i)
    if i==2: 
        force.setForceGroup(2)
        force.setReciprocalSpaceForceGroup(1)
    else: 
        force.setForceGroup(2)

platform = mm.Platform.getPlatformByName('Reference')
simulationTest = app.Simulation(pdb.topology, systemTest, integratorTest, platform)
simulationTest.context.setPositions(pdb.positions)
print "Energy: ", simulationTest.context.getState(getEnergy=True, groups=2).getPotentialEnergy()

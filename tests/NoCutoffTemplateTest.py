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

InFile = 'OneWater2.pdb'
pdb = app.PDBFile(InFile)
pdb.topology.setUnitCellDimensions((2,2,2))
forcefield = app.ForceField('tip3p.xml')

systemTest = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, 
    constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

print "this is a template TIP3P NoCutoff energy test for ", InFile
print "Note: LJ is turned off(epsilons=0) and water is set to be rigid..."

PME = systemTest.getForce(2)
N_PARTICLES = systemTest.getNumParticles()
print "#of particles: ", N_PARTICLES
for i in range(N_PARTICLES):
    charge ,sigma, epsilon = PME.getParticleParameters(i)
    PME.setParticleParameters(i, charge, sigma, 0.0)

N_EXCEPTIONS = PME.getNumExceptions()
print "#of exceptions: ", N_EXCEPTIONS
print "\n"
for i in range(N_EXCEPTIONS):
    Q1, Q2, QProd = PME.getExceptionParameters(i)[:3]
    print "Exception "+str(i)+": Q"+str(Q1)+"*Q"+str(Q2)+" = ", QProd
    print '\n'

integratorTest = mm.VerletIntegrator(0.002)
for i in range(len(systemTest.getForces())):
    force = systemTest.getForce(i)
    if i==2: 
        force.setForceGroup(1)
        force.setReciprocalSpaceForceGroup(1)
    else: 
        force.setForceGroup(2)

platform = mm.Platform.getPlatformByName('Reference')
simulationTest = app.Simulation(pdb.topology, systemTest, integratorTest, platform)
simulationTest.context.setPositions(pdb.positions)

print "Total NoCutoff Energy: ", simulationTest.context.getState(getEnergy=True, groups=2).getPotentialEnergy()
print "Forces: " 
forces = simulationTest.context.getState(getEnergy=True, getForces=True, groups=2).getForces()
for i in range(len(forces)):
    print "particle ", i, forces[i] 

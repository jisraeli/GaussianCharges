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

pdb = app.PDBFile('TwoWaters.pdb')
pdb.topology.setUnitCellDimensions((2,2,2))
forcefield = app.ForceField('tip3p.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, 
    constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)
'''
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff = 1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)
'''
PME = system.getForce(2)
PME.setUseDispersionCorrection(False)
print "Uses dispersion correciton: ", PME.getUseDispersionCorrection()
N_PARTICLES = system.getNumParticles()

for i in range(N_PARTICLES):
    charge ,sigma, epsilon = PME.getParticleParameters(i)
    PME.setParticleParameters(i, charge, sigma, 0.0)

integrator = mm.VerletIntegrator(2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)

platform = mm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

print "Energy: ", simulation.context.getState(getEnergy=True).getPotentialEnergy()
sys.exit()


'''
print('Minimizing...')
simulation.minimizeEnergy()
print simulation.context.getState(getEnergy=True).getPotentialEnergy()

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
print('Equilibrating...')
simulation.step(100)

simulation.reporters.append(app.DCDReporter('CoulombTrajectory.dcd', 2))
simulation.reporters.append(app.StateDataReporter(stdout, 10, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=1000, separator='\t'))

print('Running Production...')
simulation.step(2000)
print('Done!')
'''

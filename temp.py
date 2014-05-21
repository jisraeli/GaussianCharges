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

'''
TODO
ask Peter if script/sys is supported on openmm 6?
'''

epsilon = 8.854187817620E-12*farad/meter
COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)
CUTOFF_DIST = 1*nanometer

InFile = './examples/WaterBox.pdb'
pdb = app.PDBFile(InFile)
pdb.topology.setUnitCellDimensions((2.48,2.48,2.48))
forcefield = app.ForceField('/Users/jisraeli/Desktop/GaussianCharges/tip3pGaussian.xml')

'''
setup system
'''
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)
print "works so far"
#PME = system.getForce(2)
#N_PARTICLES = system.getNumParticles()
integrator = mm.VerletIntegrator(0.001)
integrator.setConstraintTolerance(0.00001)
thermostat = mm.AndersenThermostat(298.15*unit.kelvin, 1.0/(2.0*unit.picosecond))
barostat = mm.MonteCarloBarostat(1*unit.atmospheres, 298.15*unit.kelvin, 25)
system.addForce(thermostat)
system.addForce(barostat)
'''
add GaussianPME_DirectSpace modification with customNonBondedForce
'''
'''
forceCustomNonBonded = mm.CustomNonbondedForce("-COULOMB_CONSTANT*q1*q2*erfc(p*r)/r")
forceCustomNonBonded.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
a, b = [1.0/((0.01*nanometer)**2)]*2
p = sqrt(a * b / (a + b))
forceCustomNonBonded.addGlobalParameter("p", p)
forceCustomNonBonded.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceCustomNonBonded.addPerParticleParameter("q")
for i in range(N_PARTICLES):
    params = PME.getParticleParameters(i)[0]
    forceCustomNonBonded.addParticle([params])
'''
'''
add exceptions from TIP3P
'''
'''
for i in range(PME.getNumExceptions()):
    Q1, Q2, QProd = PME.getExceptionParameters(i)[:3]
    forceCustomNonBonded.addExclusion(Q1, Q2)
system.addForce(forceCustomNonBonded)
'''
'''
create simulation object and integrate
'''
platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(298.15*unit.kelvin)
print('Equilibrating...')
simulation.step(100000)

simulation.reporters.append(app.DCDReporter('test.dcd', 10))
simulation.reporters.append(app.StateDataReporter('test.csv', 10, step=True, 
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, 
    progress=True, remainingTime=True, density=True, speed=True, 
    totalSteps=5000, separator='\t'))
simulation.reporters.append(app.StateDataReporter(stdout, 10, step=True, 
    totalEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=5000, separator='\t'))

print('Running Production...')
simulation.step(5000)
print('Done!')




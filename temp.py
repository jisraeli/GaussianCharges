from math import*
import numpy as np
from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import*
from simtk import unit
from pylab import*
from sys import stdout

epsilon = 8.854187817620E-12*farad/meter
COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)
CUTOFF_DIST = 1*nanometer

InFile = './examples/WaterBox.pdb'
pdb = app.PDBFile(InFile)
pdb.topology.setUnitCellDimensions((2.48,2.48,2.48))
forcefield = app.ForceField('./forcefields/tip3pGaussian.xml')

'''
setup system
'''
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

integrator = mm.VerletIntegrator(0.001)
integrator.setConstraintTolerance(0.00001)
thermostat = mm.AndersenThermostat(298.15*unit.kelvin, 1.0/(2.0*unit.picosecond))
barostat = mm.MonteCarloBarostat(1*unit.atmospheres, 298.15*unit.kelvin, 25)
system.addForce(thermostat)
system.addForce(barostat)
'''
create simulation object and integrate
'''
platform = mm.Platform.getPlatformByName('OpenCL')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(298.15*unit.kelvin)
print('Equilibrating...')
simulation.step(100000)

simulation.reporters.append(app.DCDReporter('test_0.01.dcd', 100))
simulation.reporters.append(app.StateDataReporter('test_0.01.csv', 100, step=True, 
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, 
    progress=True, remainingTime=True, density=True, speed=True, 
    totalSteps=5000000, separator='\t'))
simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, 
    totalEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=5000000, separator='\t'))

print('Running Production...')
simulation.step(5000000)
print('Done!')




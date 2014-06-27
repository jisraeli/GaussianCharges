from math import*
import numpy as np
from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import*
from simtk import unit
from pylab import*
from sys import stdout
import sys

epsilon = 8.854187817620E-12*farad/meter
COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)

InFile = './WaterBoxes/WaterBox.pdb'
pdb = app.PDBFile(InFile)
pdb.topology.setUnitCellDimensions((2.48,2.48,2.48))
forcefield = app.ForceField('tip3p.xml')
'''
setup system
'''
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

PME = system.getForce(2)
N_PARTICLES = system.getNumParticles()
integrator = mm.VerletIntegrator(0.001)
integrator.setConstraintTolerance(0.00001)
thermostat = mm.AndersenThermostat(298.15*unit.kelvin, 1.0/(2.0*unit.picosecond))
barostat = mm.MonteCarloBarostat(1*unit.atmospheres, 298.15*unit.kelvin, 25)
system.addForce(thermostat)
system.addForce(barostat)
'''
add GaussianPME_DirectSpace modification with customNonBondedForce
'''
forceCustomNonBonded = mm.CustomNonbondedForce("-COULOMB_CONSTANT*q1*q2*erfc(p*r)/r; p=sqrt(w1*w2/(w1+w2))")
if PME.getNonbondedMethod() in [2, 3, 4]:
    forceCustomNonBonded.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
elif PME.getNonbondedMethod() in [1]:
    forceCustomNonBonded.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
elif PME.getNonbondedMethod() in [0]:
    forceCustomNonBonded.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)

width_H = 0.01
width_O = 0.02


def Gwidth(width):
    return 1.0/((width*nanometer)**2)

forceCustomNonBonded.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceCustomNonBonded.addPerParticleParameter("q")
forceCustomNonBonded.addPerParticleParameter("w")

# get max charge(hydrogen) and min charge(oxygen)
chargeList = []
for i in range(N_PARTICLES):
    charge = PME.getParticleParameters(i)[0]
    chargeList.append(charge)
chargeMin = min(chargeList)
chargeMax = max(chargeList)

# create widthDict dictionary so widthDict['charge'] = width of that charge
widthDict = {}
widthDict[str(chargeMax)] = Gwidth(width_H)
widthDict[str(chargeMin)] = Gwidth(width_O)

# copy charges into customForce and assign widths using dictionary
for i in range(N_PARTICLES):
    charge = PME.getParticleParameters(i)[0]
    forceCustomNonBonded.addParticle([charge, widthDict[str(charge)]])
'''
add exceptions from TIP3P
'''
for i in range(PME.getNumExceptions()):
    Q1, Q2, QProd = PME.getExceptionParameters(i)[:3]
    forceCustomNonBonded.addExclusion(Q1, Q2)
system.addForce(forceCustomNonBonded)
'''
create simulation object and integrate
'''
platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
sys.exit()
print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(298.15*unit.kelvin)
print('Equilibrating...')
simulation.step(100000)

simulation.reporters.append(app.DCDReporter('GaussianTraj.dcd', 100))
simulation.reporters.append(app.StateDataReporter('Gaussian.csv', 100, step=True, 
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, 
    progress=True, remainingTime=True, density=True, speed=True, 
    totalSteps=5000000, separator='\t'))
simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, 
    totalEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=5000000, separator='\t'))

print('Running Production...')
simulation.step(5000000)
print('Done!')
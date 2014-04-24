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

epsilon = 8.854187817620E-12*farad/meter
COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)
CUTOFF_DIST = 1*nanometer

InFile = 'TwoWaters.pdb'
pdb = app.PDBFile(InFile)
pdb.topology.setUnitCellDimensions((2,2,2))
forcefield = app.ForceField('tip3p.xml')
'''
setup system:
move all forces but PME_direct to group 1
move move PME_direct to group 2
'''
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

PME = system.getForce(2)
N_PARTICLES = system.getNumParticles()
for force in system.getForces():
    if type(force)==type(mm.NonbondedForce()): 
        force.setForceGroup(2)
        force.setReciprocalSpaceForceGroup(1)
    else: 
        force.setForceGroup(1)
integrator = mm.VerletIntegrator(0.002)
'''
add GaussianPME_DirectSpace and LJ through customNonBondedForce in group1
'''
forceCustom = mm.CustomNonbondedForce("COULOMB_CONSTANT*q1*q2*(erfc(ALPHA*r))/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
forceCustom.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
forceCustom.setForceGroup(1)

a, b = [1.0/((0.2*nanometer)**2)]*2
p = sqrt(a * b / (a + b))
ERROR_TOL = pmeCustom.getEwaldErrorTolerance()
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST

forceCustom.addGlobalParameter("ALPHA", ALPHA)
forceCustom.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceCustom.addPerParticleParameter("q")
forceCustom.addPerParticleParameter("sigma")
forceCustom.addPerParticleParameter("epsilon")
for i in range(N_PARTICLES):
    params = PME.getParticleParameters(i)
    forceCustom.addParticle(params)
'''
TODO: copy exceptions from TIP3P(temporary code below is wrong!)

for i in range(pmeTemp.getNumExceptions()):
    Q1, Q2, QProd = pmeCustom.getExceptionParameters(i)[:3]
    forceCustom.addExclusion(Q1, Q2)
'''
system.addForce(forceCustom)

platform = mm.Platform.getPlatformByName('Reference')
simulation = simulation1.Simulation1(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
'''
TODO: add dispersion correction and thermostat for equilibration
print('Equilibrating...')
simulation.step(100)
'''
simulation.reporters.append(app.DCDReporter('GaussianTraj.dcd', 5))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=1000, separator='\t'))

print('Running Production...')
simulation.step(1000)
print('Done!')



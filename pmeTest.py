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
'''
Set up system with TIP3P, forces in group 1,
and PME_Direct in group 0
'''
pdb = app.PDBFile('TwoWaters.pdb')
pdb.topology.setUnitCellDimensions((2,2,2))
forcefield = app.ForceField('tip3p.xml')
systemTest = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)
integratorTest = mm.CustomIntegrator(0.002)
integratorTest.addPerDofVariable("x1", 0)
integratorTest.addUpdateContextState();
integratorTest.addComputePerDof("v", "v+0.5*dt*f1/m")
integratorTest.addComputePerDof("x", "x+dt*v")
integratorTest.addComputePerDof("x1", "x")
integratorTest.addConstrainPositions()
integratorTest.addComputePerDof("v", "v+0.5*dt*f1/m+(x-x1)/dt")
integratorTest.addConstrainVelocities()
for i in range(len(systemTest.getForces())):
    force = systemTest.getForce(i)
    if i==2: 
        force.setReciprocalSpaceForceGroup(1)
    else: 
        force.setForceGroup(1)
'''
Create PME_direct force with leanord-Jones potential
'''
N_PARTICLES = systemTest.getNumParticles()
PME = systemTest.getForce(2)
ERROR_TOL = PME.getEwaldErrorTolerance()
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
forceTest = mm.CustomNonbondedForce("COULOMB_CONSTANT*q1*q2*(erfc(ALPHA*r))/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6)")
forceTest.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
for i in range(PME.getNumExceptions()):
    Particles = PME.getExceptionParameters(i)[:2]
    forceTest.addExclusion(Particles[0], Particles[1])
forceTest.setForceGroup(1)
forceTest.addGlobalParameter("ALPHA", ALPHA)
forceTest.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceTest.addPerParticleParameter("q")
forceTest.addPerParticleParameter("sigma")
forceTest.addPerParticleParameter("epsilon")
for i in range(N_PARTICLES):
    forceTest.addParticle(PME.getParticleParameters(i))
systemTest.addForce(forceTest)
'''
create a simulation1 object and integrate
'''
platform = mm.Platform.getPlatformByName('Reference')
simulationTest = simulation1.Simulation1(pdb.topology, systemTest, integratorTest, platform)
simulationTest.context.setPositions(pdb.positions)
print('Minimizing...')
simulationTest.minimizeEnergy()
print simulationTest.context.getState(getEnergy=True).getPotentialEnergy()
simulationTest.context.setVelocitiesToTemperature(300*unit.kelvin)
print('Equilibrating...')
simulationTest.step(100)
print('Simulating...')
simulationTest.reporters.append(app.DCDReporter('pmeTest.dcd', 2))
simulationTest.reporters.append(app.StateDataReporter(stdout, 1, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=2000, separator='\t'))
simulationTest.step(2000)


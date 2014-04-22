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
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)
integrator = mm.CustomIntegrator(0.002)
integrator.addPerDofVariable("x1", 0)
integrator.addUpdateContextState();
integrator.addComputePerDof("v", "v+0.5*dt*f1/m")
integrator.addComputePerDof("x", "x+dt*v")
integrator.addComputePerDof("x1", "x")
integrator.addConstrainPositions()
integrator.addComputePerDof("v", "v+0.5*dt*f1/m+(x-x1)/dt")
integrator.addConstrainVelocities()
for i in range(len(system.getForces())):
    force = system.getForce(i)
    if i==2: 
        force.setReciprocalSpaceForceGroup(1)
    else: 
        force.setForceGroup(1)
'''
Create gaussian force with leanord-Jones potential
TODO: UseDispersionCorrection if there's a thermostat
TODO: scale 1--4 interactions in bigger molecules
'''
N_PARTICLES = system.getNumParticles()
PME = system.getForce(2)
a, b = [1.0/((0.02*nanometer)**2)]*2
p = sqrt(a * b / (a + b))
ERROR_TOL = PME.getEwaldErrorTolerance()
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
forceGaussian = mm.CustomNonbondedForce("COULOMB_CONSTANT*q1*q2*(erf(p*r)-erf(ALPHA*r))/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6)")
forceGaussian.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
for i in range(PME.getNumExceptions()):
    Particles = PME.getExceptionParameters(i)[:2]
    forceGaussian.addExclusion(Particles[0], Particles[1])
forceGaussian.setForceGroup(1)
forceGaussian.addGlobalParameter("p", p)
forceGaussian.addGlobalParameter("ALPHA", ALPHA)
forceGaussian.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceGaussian.addPerParticleParameter("q")
forceGaussian.addPerParticleParameter("sigma")
forceGaussian.addPerParticleParameter("epsilon")
for i in range(N_PARTICLES):
    forceGaussian.addParticle(PME.getParticleParameters(i))
system.addForce(forceGaussian)
'''
create a simulation1 object and integrate
'''
platform = mm.Platform.getPlatformByName('Reference')
simulation = simulation1.Simulation1(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
print('Minimizing...')
simulation.minimizeEnergy()
print simulation.context.getState(getEnergy=True).getPotentialEnergy()
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
print('Equilibrating...')
simulation.step(100)
simulation.reporters.append(app.DCDReporter('GaussianTrajectory_0.02Width.dcd', 1))
simulation.reporters.append(app.StateDataReporter(stdout, 1, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=2000, separator='\t'))
simulation.step(2000)


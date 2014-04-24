from math import*
import numpy as np
from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import*
from pylab import*

epsilon = 8.854187817620E-12*farad/meter
COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)
CUTOFF_DIST = 1.0*nanometers
'''
Get particle parameters
'''
InFile = 'TwoWaters.pdb'
pdb = app.PDBFile(InFile)
pdb.topology.setUnitCellDimensions((2,2,2))
forcefield = app.ForceField('tip3p.xml')

systemTest = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)
PME = systemTest.getForce(2)
N_PARTICLES = systemTest.getNumParticles()
'''
Setup systems and integrators
'''
system = mm.System()
systemGaussian = mm.System()
MASS = 1.0 * atomic_mass_unit   
for i in range(N_PARTICLES):
    system.addParticle(MASS)
    systemGaussian.addParticle(MASS)
integrator = mm.LangevinIntegrator(300, 1.0, 0.002)
integratorGaussian = mm.LangevinIntegrator(300, 1.0, 0.002)
'''
Set up NonBonded force with PME:
direct space in group 1
reciprocal space in group 2
'''
force = mm.NonbondedForce()
force.setNonbondedMethod(mm.NonbondedForce.PME)
force.setCutoffDistance(CUTOFF_DIST)
force.setForceGroup(1)
force.setReciprocalSpaceForceGroup(2)
for i in range(N_PARTICLES):
    charge ,sigma, epsilon = PME.getParticleParameters(i)
    force.addParticle(charge, sigma, 0.0)
system.addForce(force)
'''
Setup customForce with direct space part of PME and LJ
'''
a, b = [1.0/((0.4*nanometer)**2)]*2
p = sqrt(a * b / (a + b))
ERROR_TOL = force.getEwaldErrorTolerance()
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
forceGaussian = mm.CustomNonbondedForce("COULOMB_CONSTANT*q1*q2*erfc(ALPHA*r)/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
forceGaussian.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
forceGaussian.setCutoffDistance(CUTOFF_DIST)
forceGaussian.setForceGroup(1)
forceGaussian.addGlobalParameter("ALPHA", ALPHA)
forceGaussian.addGlobalParameter("p", p)
forceGaussian.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceGaussian.addPerParticleParameter("q")
forceGaussian.addPerParticleParameter("sigma")
forceGaussian.addPerParticleParameter("epsilon")
for i in range(N_PARTICLES):
    charge ,sigma, epsilon = PME.getParticleParameters(i)
    forceGaussian.addParticle([charge, sigma, 0.0])
systemGaussian.addForce(forceGaussian)
'''
setup contexts and print out direct space forces
'''
context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))

contextGaussian = mm.Context(systemGaussian, integratorGaussian, mm.Platform.getPlatformByName('Reference'))
context.setPositions(pdb.positions)
contextGaussian.setPositions(pdb.positions)
state = context.getState(getEnergy=True, getPositions=True, getForces=True, groups=2)
stateGaussian = contextGaussian.getState(getEnergy=True, getPositions=True, getForces=True, groups=2)
print "NonBondedForce direct space PME forces: \n"
forces = state.getForces()
for i in range(len(forces)):
    print "particle ", i, forces[i]
print "CustomNonbondedForce direct space PME forces: \n"
forcesGaussian = stateGaussian.getForces()
for i in range(len(forces)):
    print "particle ", i, forcesGaussian[i]






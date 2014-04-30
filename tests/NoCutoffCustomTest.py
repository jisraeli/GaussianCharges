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

InFile = 'OneWater2.pdb'
pdb = app.PDBFile(InFile)
pdb.topology.setUnitCellDimensions((2,2,2))
forcefield = app.ForceField('tip3p.xml')

systemTest = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, 
    constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

print "this is a custom TIP3P PME energy test for ", InFile
print "Note: LJ is turned off(epsilons=0) and water is set to be rigid..."
print "Setting up system..."
PME = systemTest.getForce(2)
N_PARTICLES = systemTest.getNumParticles()
for i in range(N_PARTICLES):
    charge ,sigma, epsilon = PME.getParticleParameters(i)
    PME.setParticleParameters(i, charge, sigma, 0.0)
print "Isolating default PME..."
integratorTest = mm.VerletIntegrator(0.002)
for i in range(len(systemTest.getForces())):
    force = systemTest.getForce(i)
    if i==2: 
        force.setForceGroup(2)
        force.setReciprocalSpaceForceGroup(2)
    else: 
        force.setForceGroup(2)
print "Creating Direct space PME through customForce..."
N_PARTICLES = systemTest.getNumParticles()
PME = systemTest.getForce(2)
ERROR_TOL = PME.getEwaldErrorTolerance()
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
forceTest = mm.CustomNonbondedForce("COULOMB_CONSTANT*q1*q2*erfc(ALPHA*r)/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
forceTest.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
print "Copying TIP3P's exceptions: \n"
for i in range(PME.getNumExceptions()):
    Q1, Q2, QProd = PME.getExceptionParameters(i)[:3]
    print "Exception "+str(i)+": Q"+str(Q1)+"*Q"+str(Q2)+" = ", QProd
    forceTest.addExclusion(Q1, Q2)
forceTest.setForceGroup(1)
forceTest.addGlobalParameter("ALPHA", ALPHA)
forceTest.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceTest.addPerParticleParameter("q")
forceTest.addPerParticleParameter("sigma")
forceTest.addPerParticleParameter("epsilon")
ESELF = 0
for i in range(N_PARTICLES):
    forceTest.addParticle(PME.getParticleParameters(i))
    ESELF += PME.getParticleParameters(i)[0].value_in_unit(elementary_charge)**2
ESELF *= -ALPHA.value_in_unit(ALPHA.unit)*COULOMB_CONSTANT/sqrt(pi)
systemTest.addForce(forceTest)
'''
create a simulation1 object and integrate
'''
platform = mm.Platform.getPlatformByName('Reference')
simulationTest = app.Simulation(pdb.topology, systemTest, integratorTest, platform)
simulationTest.context.setPositions(pdb.positions)

print "\nDirect Space Custom PME Energy: ", simulationTest.context.getState(getEnergy=True, groups=2).getPotentialEnergy()
print "Custom PME Self Energy: ", ESELF*kilojoule_per_mole
print "Forces: " 
forces = simulationTest.context.getState(getEnergy=True, getForces=True, groups=2).getForces()
for i in range(len(forces)):
    print "particle ", i, forces[i] 


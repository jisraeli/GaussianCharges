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

InFile = 'WaterBox.pdb'
pdb = app.PDBFile(InFile)
pdb.topology.setUnitCellDimensions((2,2,2))
forcefield = app.ForceField('tip3p.xml')
'''
setup systemTemp, systemCustomBond, systemCustomNonBonded
'''
systemTemp = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

systemCustomNonBonded = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

pmeTemp = systemTemp.getForce(2)
pmeCustomNonBonded = systemCustomNonBonded.getForce(2)
print "Uses dispersion correction: ", pmeTemp.getUseDispersionCorrection()
'''
make list of TIP3P's exception pairs
'''
ExceptionPairs = []
for i in range(pmeTemp.getNumExceptions()):
    Q1, Q2, QProd = pmeTemp.getExceptionParameters(i)[:3]
    ExceptionPairs.append([Q1, Q2])
'''
Turn Off LJ in every system
move systemCustomBond forces to group 2
move systemCustomNonBonded forces to group 2
move systemTemp forces to group 2, but leave PME_direct in group 1
'''
N_PARTICLES = systemTemp.getNumParticles()

for i in range(len(systemTemp.getForces())):
    force = systemTemp.getForce(i)
    if i==2: 
        force.setForceGroup(1)
        force.setReciprocalSpaceForceGroup(1)
    else: 
        force.setForceGroup(1)

for i in range(len(systemCustomNonBonded.getForces())):
    force = systemCustomNonBonded.getForce(i)
    if i==2: 
        force.setForceGroup(2)
        force.setReciprocalSpaceForceGroup(1)
    else: 
        force.setForceGroup(1)
integratorTemp = mm.VerletIntegrator(0.002)
integratorCustomNonBonded = mm.VerletIntegrator(0.002)
'''
add PME_direct to systemCustomBond w CustomBondForce
add TIP3P's exceptions by setting those bonds with
    charge_product=0
'''
ERROR_TOL = pmeTemp.getEwaldErrorTolerance()
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
forceCustomBond = mm.CustomBondForce("-COULOMB_CONSTANT*Q*erf(ALPHA*r)/r")
forceCustomBond.setForceGroup(1)
forceCustomBond.addGlobalParameter("ALPHA", ALPHA)
forceCustomBond.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceCustomBond.addPerBondParameter("Q")
for i in range(N_PARTICLES):
    for j in range(i+1, N_PARTICLES):
        if [i, j] in ExceptionPairs:
            Q1 ,sigma1, epsilon1 = pmeTemp.getParticleParameters(i)
            Q2 ,sigma2, epsilon2 = pmeTemp.getParticleParameters(j)
            forceCustomBond.addBond(i, j, [Q1*Q2])
systemCustomNonBonded.addForce(forceCustomBond)
'''
add PME_direct to systemCustomNonBonded w CustomNonbondedForce
add exception pairs from TIP3P as exclusions
'''
forceCustomNonBonded = mm.CustomNonbondedForce("COULOMB_CONSTANT*q1*q2*erfc(ALPHA*r)/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6);sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
forceCustomNonBonded.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
forceCustomNonBonded.setUseLongRangeCorrection(True)
print "customNonBonded uses long range correction: ", forceCustomNonBonded.getUseLongRangeCorrection()
forceCustomNonBonded.setForceGroup(1)
forceCustomNonBonded.addGlobalParameter("ALPHA", ALPHA)
forceCustomNonBonded.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceCustomNonBonded.addPerParticleParameter("q")
forceCustomNonBonded.addPerParticleParameter("sigma")
forceCustomNonBonded.addPerParticleParameter("epsilon")
for i in range(N_PARTICLES):
    charge ,sigma, epsilon = pmeTemp.getParticleParameters(i)
    forceCustomNonBonded.addParticle([charge, sigma, epsilon])

for i in range(pmeTemp.getNumExceptions()):
    Q1, Q2, QProd = pmeTemp.getExceptionParameters(i)[:3]
    forceCustomNonBonded.addExclusion(Q1, Q2)
systemCustomNonBonded.addForce(forceCustomNonBonded)
'''
create simulation objects and print forces
'''
platformTemp = mm.Platform.getPlatformByName('Reference')
platformCustomNonBonded = mm.Platform.getPlatformByName('Reference')
simulationTemp = app.Simulation(pdb.topology, systemTemp, integratorTemp, platformTemp)
simulationTemp.context.setPositions(pdb.positions)
simulationCustomNonBonded = app.Simulation(pdb.topology, systemCustomNonBonded, integratorCustomNonBonded, platformCustomNonBonded)
simulationCustomNonBonded.context.setPositions(pdb.positions)

def RelativeForceError(force, forceCustom):
    force = np.array(force)
    forceCustom = np.array(forceCustom)
    diff = force - forceCustom
    RelativeDiff = diff/np.linalg.norm(force)
    return np.linalg.norm(RelativeDiff)

print "\nNonBondedForce Total Initial Energy: " , simulationTemp.context.getState(getEnergy=True, groups=2).getPotentialEnergy() 
forces = simulationTemp.context.getState(getEnergy=True, getForces=True, groups=2).getForces()

print "\nCustomNonBondedForce Total Initial Energy: " , simulationCustomNonBonded.context.getState(getEnergy=True, groups=2).getPotentialEnergy() 
forcesCustom = simulationCustomNonBonded.context.getState(getEnergy=True, getForces=True, groups=2).getForces()
RelErrList = []
for i in range(len(forces)):
    RelErr = RelativeForceError(forces[i], forcesCustom[i])
    RelErrList.append(RelErr)
RelErrVector = np.array(RelErrList)
print "\naverage relative error in initial forces: ", np.average(RelErrVector)
print "\nstandard deviation of erros in forces: ", np.std(RelErrVector)
    

'''
OUTPUT

NonBondedForce Total Initial Energy:  3079.44138935 kJ/mol

CustomNonBondedForce Total Initial Energy:  3096.39468381 kJ/mol

average relative error in initial forces:  4.31718844703e-09

standard deviation of erros in forces:  2.00889151415e-09

'''

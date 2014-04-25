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
setup systemTemp, systemCustomBond, systemCustomNonBonded
'''
systemTemp = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

systemCustomBond = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

systemCustomNonBonded = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

pmeTemp = systemTemp.getForce(2)
pmeCustomBond = systemCustomBond.getForce(2)
pmeCustomNonBonded = systemCustomNonBonded.getForce(2)
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

for i in range(N_PARTICLES):
    charge ,sigma, epsilon = pmeTemp.getParticleParameters(i)
    pmeTemp.setParticleParameters(i, charge, sigma, 0.0)
    pmeCustomBond.setParticleParameters(i, charge, sigma, 0.0)
    pmeCustomNonBonded.setParticleParameters(i, charge, sigma, 0.0)

for i in range(len(systemTemp.getForces())):
    force = systemTemp.getForce(i)
    if i==2: 
        force.setForceGroup(1)
        force.setReciprocalSpaceForceGroup(2)
    else: 
        force.setForceGroup(2)

for i in range(len(systemCustomBond.getForces())):
    force = systemCustomBond.getForce(i)
    if i==2: 
        force.setForceGroup(2)
        force.setReciprocalSpaceForceGroup(2)
    else: 
        force.setForceGroup(2)

for i in range(len(systemCustomNonBonded.getForces())):
    force = systemCustomNonBonded.getForce(i)
    if i==2: 
        force.setForceGroup(2)
        force.setReciprocalSpaceForceGroup(2)
    else: 
        force.setForceGroup(2)
integratorTemp = mm.VerletIntegrator(0.002)
integratorCustomBond = mm.VerletIntegrator(0.002)
integratorCustomNonBonded = mm.VerletIntegrator(0.002)
'''
add PME_direct to systemCustomBond w CustomBondForce
add TIP3P's exceptions by setting those bonds with
    charge_product=0
'''
ERROR_TOL = pmeTemp.getEwaldErrorTolerance()
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
forceCustomBond = mm.CustomBondForce("COULOMB_CONSTANT*Q*erfc(ALPHA*r)/r")
forceCustomBond.setForceGroup(1)
forceCustomBond.addGlobalParameter("ALPHA", ALPHA)
forceCustomBond.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceCustomBond.addPerBondParameter("Q")
for i in range(N_PARTICLES):
    for j in range(i+1, N_PARTICLES):
        if [i, j] in ExceptionPairs:
            forceCustomBond.addBond(i, j, [0*elementary_charge**2])
        else:
            Q1 ,sigma1, epsilon1 = pmeTemp.getParticleParameters(i)
            Q2 ,sigma2, epsilon2 = pmeTemp.getParticleParameters(j)
            forceCustomBond.addBond(i, j, [Q1*Q2])
systemCustomBond.addForce(forceCustomBond)
'''
add PME_direct to systemCustomNonBonded w CustomNonbondedForce
add exception pairs from TIP3P as exclusions
'''
forceCustomNonBonded = mm.CustomNonbondedForce("COULOMB_CONSTANT*q1*q2*erfc(ALPHA*r)/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6);sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
forceCustomNonBonded.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
forceCustomNonBonded.setForceGroup(1)
forceCustomNonBonded.addGlobalParameter("ALPHA", ALPHA)
forceCustomNonBonded.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceCustomNonBonded.addPerParticleParameter("q")
forceCustomNonBonded.addPerParticleParameter("sigma")
forceCustomNonBonded.addPerParticleParameter("epsilon")
for i in range(N_PARTICLES):
    charge ,sigma, epsilon = pmeTemp.getParticleParameters(i)
    forceCustomNonBonded.addParticle([charge, sigma, 0.0])

for i in range(pmeTemp.getNumExceptions()):
    Q1, Q2, QProd = pmeTemp.getExceptionParameters(i)[:3]
    forceCustomNonBonded.addExclusion(Q1, Q2)
systemCustomNonBonded.addForce(forceCustomNonBonded)
'''
create simulation objects and print forces
'''
platformTemp = mm.Platform.getPlatformByName('Reference')
platformCustomBond = mm.Platform.getPlatformByName('Reference')
platformCustomNonBonded = mm.Platform.getPlatformByName('Reference')
simulationTemp = app.Simulation(pdb.topology, systemTemp, integratorTemp, platformTemp)
simulationTemp.context.setPositions(pdb.positions)
simulationCustomBond = app.Simulation(pdb.topology, systemCustomBond, integratorCustomBond, platformCustomBond)
simulationCustomBond.context.setPositions(pdb.positions)
simulationCustomNonBonded = app.Simulation(pdb.topology, systemCustomNonBonded, integratorCustomNonBonded, platformCustomNonBonded)
simulationCustomNonBonded.context.setPositions(pdb.positions)

print "\nNonBondedForce Direct Space PME Forces: \n" 
forces = simulationTemp.context.getState(getEnergy=True, getForces=True, groups=2).getForces()
for i in range(len(forces)):
    print "particle ", i, forces[i]

print "\nCustomNonBondedForce Direct Space PME Forces: \n" 
forces = simulationCustomNonBonded.context.getState(getEnergy=True, getForces=True, groups=2).getForces()
for i in range(len(forces)):
    print "particle ", i, forces[i]

print "\nCustomBondForce Direct Space PME Forces: \n" 
forces = simulationCustomBond.context.getState(getEnergy=True, getForces=True, groups=2).getForces()
for i in range(len(forces)):
    print "particle ", i, forces[i]

'''
OUTPUT

NonBondedForce Direct Space PME Forces: 

particle  0 (-22.709305459133372, 26.97954925345011, -30.759053554131224) kJ/(nm mol)
particle  1 (-10.310511666994763, -6.259688314276094, 23.51882350862853) kJ/(nm mol)
particle  2 (3.8780704284737837, -10.966394591058279, 11.980985667127623) kJ/(nm mol)
particle  3 (54.4543816221644, -3.7010536100475804, 28.940374738738324) kJ/(nm mol)
particle  4 (-5.216419237506585, -9.99127057596625, -24.067745799572055) kJ/(nm mol)
particle  5 (-20.096215687003458, 3.938857837898105, -9.613384560791182) kJ/(nm mol)

CustomNonBondedForce Direct Space PME Forces: 

particle  0 (-69.73759608537662, 19.356122167791845, 27.06323194271893) kJ/(nm mol)
particle  1 (15.156999250665127, -3.2352719838242994, -3.9047738860232766) kJ/(nm mol)
particle  2 (25.438850003478343, -6.367383791144251, -18.417702413340155) kJ/(nm mol)
particle  3 (69.65131524300209, -39.291716675119716, -35.23000053289787) kJ/(nm mol)
particle  4 (-15.335964738338324, 7.137110106040673, 7.809693598358306) kJ/(nm mol)
particle  5 (-25.173603673430613, 22.40114017625575, 22.679551291184058) kJ/(nm mol)

CustomBondForce Direct Space PME Forces: 

particle  0 (-69.73759608537662, 19.356122167791845, 27.06323194271893) kJ/(nm mol)
particle  1 (15.156999250665127, -3.2352719838242994, -3.9047738860232766) kJ/(nm mol)
particle  2 (25.438850003478343, -6.367383791144251, -18.417702413340155) kJ/(nm mol)
particle  3 (69.65131524300209, -39.291716675119716, -35.23000053289787) kJ/(nm mol)
particle  4 (-15.335964738338324, 7.137110106040673, 7.809693598358306) kJ/(nm mol)
particle  5 (-25.173603673430613, 22.40114017625575, 22.679551291184058) kJ/(nm mol)
'''

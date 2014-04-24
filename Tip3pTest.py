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
setup systemTemp and systemCustom:
move systemCustom forces to group 2
move systemTemp forces to group 2, but leave PME_direct in group 1
'''
systemTemp = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

systemCustom = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

pmeTemp = systemTemp.getForce(2)
pmeCustom = systemCustom.getForce(2)
N_PARTICLES = systemTemp.getNumParticles()

for i in range(N_PARTICLES):
    charge ,sigma, epsilon = pmeTemp.getParticleParameters(i)
    pmeTemp.setParticleParameters(i, charge, sigma, 0.0)
    pmeCustom.setParticleParameters(i, charge, sigma, 0.0)

for i in range(len(systemTemp.getForces())):
    force = systemTemp.getForce(i)
    if i==2: 
        force.setForceGroup(1)
        force.setReciprocalSpaceForceGroup(2)
    else: 
        force.setForceGroup(2)

for i in range(len(systemCustom.getForces())):
    force = systemCustom.getForce(i)
    if i==2: 
        force.setForceGroup(2)
        force.setReciprocalSpaceForceGroup(2)
    else: 
        force.setForceGroup(2)
integratorTemp = mm.VerletIntegrator(0.002)
integratorCustom = mm.VerletIntegrator(0.002)
'''
introduce PME direct space through customNonBondedForce into systemCustom
'''
ERROR_TOL = pmeCustom.getEwaldErrorTolerance()
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
forceCustom = mm.CustomNonbondedForce("COULOMB_CONSTANT*q1*q2*erfc(ALPHA*r)/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
forceCustom.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
'''
copy exceptions from TIP3P into forceCustom as exclusions(exceptions not available)
'''
for i in range(pmeTemp.getNumExceptions()):
    Q1, Q2, QProd = pmeCustom.getExceptionParameters(i)[:3]
    forceCustom.addExclusion(Q1, Q2)
forceCustom.setForceGroup(1)
forceCustom.addGlobalParameter("ALPHA", ALPHA)
forceCustom.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceCustom.addPerParticleParameter("q")
forceCustom.addPerParticleParameter("sigma")
forceCustom.addPerParticleParameter("epsilon")
for i in range(N_PARTICLES):
    charge ,sigma, epsilon = pmeCustom.getParticleParameters(i)
    forceCustom.addParticle([charge, sigma, 0.0])
systemCustom.addForce(forceCustom)

platformTemp = mm.Platform.getPlatformByName('Reference')
platformCustom = mm.Platform.getPlatformByName('Reference')
simulationTemp = app.Simulation(pdb.topology, systemTemp, integratorTemp, platformTemp)
simulationTemp.context.setPositions(pdb.positions)
simulationCustom = app.Simulation(pdb.topology, systemCustom, integratorCustom, platformCustom)
simulationCustom.context.setPositions(pdb.positions)
'''
Serialize systems and print out direct space forces
'''
wfile = open("CustomSystem.txt", 'w')
wfile.write(mm.XmlSerializer.serialize(systemCustom))
wfile.close()
wfile = open("TemplateSystem.txt", 'w')
wfile.write(mm.XmlSerializer.serialize(systemTemp))
wfile.close()

print "\nNonBondedForce Direct Space PME Forces: " 
forces = simulationTemp.context.getState(getEnergy=True, getForces=True, groups=2).getForces()
for i in range(len(forces)):
    print "particle ", i, forces[i]

print "\nCustomNonBondedForce Direct Space PME Forces: " 
forces = simulationCustom.context.getState(getEnergy=True, getForces=True, groups=2).getForces()
for i in range(len(forces)):
    print "particle ", i, forces[i]



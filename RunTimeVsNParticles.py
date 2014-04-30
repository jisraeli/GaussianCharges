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

def CreateGaussianSimulation(InFile, GaussianWidth, Gaussian):
    '''
    setup system:
    move all forces but PME_direct to group 1
    move PME_direct to group 2
    '''
    pdb = app.PDBFile(InFile)
    pdb.topology.setUnitCellDimensions((2,2,2))
    forcefield = app.ForceField('tip3p.xml')

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
        nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
        ewaldErrorTolerance=0.0005)
    if not Gaussian:
        integrator = mm.VerletIntegrator(0.002)
        platform = mm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        return simulation
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
    forceCustomNonBonded = mm.CustomNonbondedForce("COULOMB_CONSTANT*q1*q2*(erf(p*r)-erf(ALPHA*r))/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
    forceCustomNonBonded.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    forceCustomNonBonded.setUseLongRangeCorrection(True)
    forceCustomNonBonded.setForceGroup(1)
    a, b = [1.0/((GaussianWidth)**2)]*2
    p = sqrt(a * b / (a + b))
    ERROR_TOL = PME.getEwaldErrorTolerance()
    ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
    forceCustomNonBonded.addGlobalParameter("ALPHA", ALPHA)
    forceCustomNonBonded.addGlobalParameter("p", p)
    forceCustomNonBonded.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
    forceCustomNonBonded.addPerParticleParameter("q")
    forceCustomNonBonded.addPerParticleParameter("sigma")
    forceCustomNonBonded.addPerParticleParameter("epsilon")
    for i in range(N_PARTICLES):
        params = PME.getParticleParameters(i)
        forceCustomNonBonded.addParticle(params)
    '''
    add exceptions for direct space and LJ
    '''
    for i in range(PME.getNumExceptions()):
        Q1, Q2, QProd = PME.getExceptionParameters(i)[:3]
        forceCustomNonBonded.addExclusion(Q1, Q2)
    system.addForce(forceCustomNonBonded)
    '''
    add reciprocal space exceptions through customBondForce in group1
    '''
    forceCustomBond = mm.CustomBondForce("-COULOMB_CONSTANT*Q*erf(ALPHA*r)/r")
    forceCustomBond.setForceGroup(1)
    forceCustomBond.addGlobalParameter("ALPHA", ALPHA)
    forceCustomBond.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
    forceCustomBond.addPerBondParameter("Q")
    '''
    get list of particle pairs in exceptions
    add bonds for those pairs
    '''
    ExceptionPairs = []
    for i in range(PME.getNumExceptions()):
        Q1, Q2, QProd = PME.getExceptionParameters(i)[:3]
        ExceptionPairs.append([Q1, Q2])

    for i in range(N_PARTICLES):
        for j in range(i+1, N_PARTICLES):
            if [i, j] in ExceptionPairs:
                Q1 ,sigma1, epsilon1 = PME.getParticleParameters(i)
                Q2 ,sigma2, epsilon2 = PME.getParticleParameters(j)
                forceCustomBond.addBond(i, j, [Q1*Q2])
    system.addForce(forceCustomBond)
    '''
    create simulation1 object and integrate
    '''
    platform = mm.Platform.getPlatformByName('Reference')
    simulation = simulation1.Simulation1(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    return simulation


def Value(quantity):
    return quantity.value_in_unit(quantity.unit)


BoxSizeList = np.linspace(1.1, 2.3, 13)
GaussianWidth = 0.2
ParticleNumberList = []
GaussianRunTimeList = []
ReferenceRunTimeList = []
for size in BoxSizeList:
    InFile = './examples/' + (str(size)+'_')*3 + 'WaterBox.pdb'
    print InFile

plot(Widths, EnergyDiffList, label="Energy Difference")
plot(Widths, AvgForceDiffList, label="Median Per-Particle Force Difference")
legend(loc=4)
xlabel("Gaussian Width (nm)")
ylabel("Relative Difference")
title("Point charges vs Gaussian Charges")
show()




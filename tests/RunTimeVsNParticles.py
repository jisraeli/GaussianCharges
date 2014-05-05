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
import time

epsilon = 8.854187817620E-12*farad/meter
COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)
CUTOFF_DIST = 1*nanometer

def CreateGaussianSimulation(InFile, size, GaussianWidth, Gaussian):
    '''
    setup system:
    move all forces but PME_direct to group 1
    move PME_direct to group 2
    '''
    pdb = app.PDBFile(InFile)
    pdb.topology.setUnitCellDimensions((size,size,size))
    forcefield = app.ForceField('tip3p.xml')

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
        nonbondedCutoff=1.0*nanometers, constraints=app.HBonds, rigidWater=True, 
        ewaldErrorTolerance=0.0005)
    if not Gaussian:
        PME = system.getForce(2)
        N_PARTICLES = system.getNumParticles()
        integrator = mm.VerletIntegrator(0.001)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(298.15*unit.kelvin, 1.0/(2.0*unit.picosecond))
        barostat = mm.MonteCarloBarostat(1*unit.atmospheres, 298.15*unit.kelvin, 25)
        system.addForce(thermostat)
        system.addForce(barostat)
        platform = mm.Platform.getPlatformByName('CPU')
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        return [simulation, N_PARTICLES]
    PME = system.getForce(2)
    N_PARTICLES = system.getNumParticles()
    for force in system.getForces():
        if type(force)==type(mm.NonbondedForce()): 
            force.setForceGroup(2)
            force.setReciprocalSpaceForceGroup(1)
        else: 
            force.setForceGroup(1)
    integrator = mm.VerletIntegrator(0.001)
    integrator.setConstraintTolerance(0.00001)
    thermostat = mm.AndersenThermostat(298.15*unit.kelvin, 1.0/(2.0*unit.picosecond))
    thermostat.setForceGroup(1)
    barostat = mm.MonteCarloBarostat(1*unit.atmospheres, 298.15*unit.kelvin, 25)
    barostat.setForceGroup(1)
    system.addForce(thermostat)
    system.addForce(barostat)
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
    platform = mm.Platform.getPlatformByName('CPU')
    simulation = simulation1.Simulation1(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    return [simulation, N_PARTICLES]

BoxSizeList = np.linspace(2.1, 2.5, 5)
GaussianWidth = 0.2
ParticleNumberList = []
GaussianRunTimeList = []
ReferenceRunTimeList = []
RatioList = []
STEPS = 100
print('\n#Particles <GaussianTime> <ReferenceTime> <Ratio>')
for size in BoxSizeList:
    InFile = '../examples/' + (str(size)+'_')*3 + 'WaterBox.pdb'
    [simulation, N_PARTICLES] = CreateGaussianSimulation(InFile, size, GaussianWidth, Gaussian=True)
    ParticleNumberList.append([N_PARTICLES])
    total = 0
    for i in range(3):
        GaussianStartTime = time.time()
        simulation.step(STEPS)
        GaussianEndTime = time.time()
        total += GaussianEndTime - GaussianStartTime
    GaussianRunTime = total/3
    [ReferenceSimulation, _] = CreateGaussianSimulation(InFile, size, GaussianWidth, Gaussian=False)
    total = 0
    for i in range(3):
        ReferenceStartTime = time.time()
        ReferenceSimulation.step(STEPS)
        ReferenceEndTime = time.time()
        total += ReferenceEndTime - ReferenceStartTime
    ReferenceRunTime = total/3
    GaussianRunTimeList.append([GaussianRunTime])
    ReferenceRunTimeList.append([ReferenceRunTime])
    RatioList.append([GaussianRunTime/ReferenceRunTime])
    column1 = np.asarray(ParticleNumberList)
    column2 = np.asarray(GaussianRunTimeList)
    column3 = np.asarray(ReferenceRunTimeList)
    column4 = np.asarray(RatioList)
    data = np.hstack((column1, column2, column3, column4))
    np.savetxt('RunTimeVsNParticles_CPU.txt', data, header="Steps <GaussianRunTime> <ReferenceRunTime> Ratio")
    print N_PARTICLES, GaussianRunTime, ReferenceRunTime, GaussianRunTime/ReferenceRunTime

plot(ParticleNumberList, RatioList)
xlabel("# of Particles")
ylabel("Gaussian Run Time / Reference PME RunTime")
title("100 Integrator Steps RunTime Comparison (CPU)")
show()
'''
OUTPUT for STEPS=100:

#Particles <GaussianTime> <ReferenceTime> <Ratio>
6 0.239636023839 0.235470612844 1.01768972759
12 0.258547465007 0.23966105779 1.078804656
27 0.278131961823 0.253934780757 1.09528895961
42 0.334018627803 0.248775641123 1.34265005325
117 0.576965411504 0.328543980916 1.75612838773
198 1.33016228676 0.524753093719 2.53483457778
267 1.83385872841 0.574059724808 3.19454344062
369 3.8232913812 0.951488018036 4.01822336039
483 5.77818369865 1.41246700287 4.09084508659
594 8.95749791463 1.88364664714 4.755402468
735 11.6400163174 2.37498164177 4.90109738646
867 15.4102863471 3.12877090772 4.92534826025
1041 22.2165226142 4.27143033346 5.20119043969
1209 29.1707986196 5.57625468572 5.2312529222

'''



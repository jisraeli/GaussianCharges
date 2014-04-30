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
        PME = system.getForce(2)
        N_PARTICLES = system.getNumParticles()
        integrator = mm.VerletIntegrator(0.002)
        platform = mm.Platform.getPlatformByName('Reference')
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
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    return [simulation, N_PARTICLES]


def Value(quantity):
    return quantity.value_in_unit(quantity.unit)


BoxSizeList = np.linspace(1.2, 2.5, 14)
GaussianWidth = 0.2
ParticleNumberList = []
GaussianRunTimeList = []
ReferenceRunTimeList = []
RatioList = []
STEPS = 100
print('\n#Particles <GaussianTime> <ReferenceTime> <Ratio>')
for size in BoxSizeList:
    InFile = '../examples/' + (str(size)+'_')*3 + 'WaterBox.pdb'
    [simulation, N_PARTICLES] = CreateGaussianSimulation(InFile, GaussianWidth, Gaussian=True)
    ParticleNumberList.append([N_PARTICLES])
    total = 0
    for i in range(3):
        GaussianStartTime = time.time()
        simulation.step(STEPS)
        GaussianEndTime = time.time()
        total += GaussianEndTime - GaussianStartTime
    GaussianRunTime = total/3
    [ReferenceSimulation, _] = CreateGaussianSimulation(InFile, GaussianWidth, Gaussian=False)
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
    np.savetxt('RunTimeVsNParticles.txt', data, header="Steps <GaussianRunTime> <ReferenceRunTime> Ratio")
    print N_PARTICLES, GaussianRunTime, ReferenceRunTime, GaussianRunTime/ReferenceRunTime

#plot(ParticleNumberList, GaussianRunTimeList, label="Gaussian")
#plot(ParticleNumberList, ReferenceRunTimeList, label="Reference")
plot(ParticleNumberList, RatioList)
#legend(loc=1)
xlabel("# of Particles")
ylabel("Gaussian Run Time / Reference PME RunTime")
title("100 Integrator Steps RunTime Comparison")
show()
'''
OUTPUT for STEPS=1000:
#Particles GaussianTime ReferenceTime
6 2.94660186768 2.62219500542
12 2.7828309536 2.92966508865
27 3.72357606888 2.62539601326
42 4.46127104759 3.44416904449
117 6.65421199799 3.33271002769
198 13.0288012028 5.79758000374
267 23.6940391064 6.42616200447
369 36.0542321205 11.3874499798
483 71.6805109978 14.8495919704
594 92.2409329414 18.2430429459
735 165.435776949 28.1915528774
867 276.529452085 45.3385272026

OUTPUT for STEPS=100:
#Particles GaussianTime ReferenceTime Ratio
6 0.248203992844 0.237093925476 1.04685935055
12 0.250403881073 0.259467124939 0.96506977958
27 0.265575170517 0.243100166321 1.09245161999
42 0.292719125748 0.247663021088 1.18192503856
117 1.08535599709 0.459707975388 2.36096838688
198 1.64134907722 0.650513887405 2.52315762815
267 2.71217298508 0.871045827866 3.11369723419
369 4.91774892807 1.42550992966 3.44981737816
483 6.71047115326 1.40629911423 4.7717239422
594 8.40508699417 1.84828805923 4.54749840112
735 13.4106059074 3.93785309792 3.40556277088
867 18.7825930119 3.92309784889 4.78769424962
1041 24.3472688198 4.30676603317 5.65326015676
1209 32.9583649635 6.09436011314 5.40801074299

Avg OUTPUT for STEPS=100:
TODO: further testing if ratio converges
'''



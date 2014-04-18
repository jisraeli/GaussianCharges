from math import*
import numpy as np
from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import*
from pylab import*

epsilon = 8.854187817620E-12*farad/meter

COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)

CUTOFF_DIST = 1*nanometer

######### Set system&Integrator #######
system = mm.System()
#BOX_LENGTH = 20
#system.setDefaultPeriodicBoxVectors((BOX_LENGTH, 0, 0), (0, BOX_LENGTH, 0), (0, 0, BOX_LENGTH))

systemGaussian = mm.System()
#systemGaussian.setDefaultPeriodicBoxVectors((BOX_LENGTH, 0, 0), (0, BOX_LENGTH, 0), (0, 0, BOX_LENGTH))


N_PARTICLES = 2
MASS = 1.0 * atomic_mass_unit   
CHARGE = [1.0, -1.0] * elementary_charge
for i in range(N_PARTICLES):
    system.addParticle(MASS)
    systemGaussian.addParticle(MASS)
integrator = mm.LangevinIntegrator(300, 1.0, 0.002)
integratorGaussian = mm.LangevinIntegrator(300, 1.0, 0.002)
######### PME #########################
force = mm.NonbondedForce()
force.setNonbondedMethod(mm.NonbondedForce.PME)
force.setCutoffDistance(CUTOFF_DIST)
force.setReciprocalSpaceForceGroup(1)
for i in range(N_PARTICLES):
    force.addParticle(CHARGE[i], 1.0, 0.0)
system.addForce(force)
######### Gaussian Potential ##########
a, b = [1.0/((0.4*nanometer)**2)]*2
p = sqrt(a * b / (a + b))
# force.setEwaldErrorTolerance(0.1)
ERROR_TOL = force.getEwaldErrorTolerance()
print "ERROR_TOL: ", ERROR_TOL
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
# forceGaussian = mm.CustomNonbondedForce("COULOMB_CONSTANT * q1 * q2 * erf(p*r) / r")
forceGaussian = mm.CustomNonbondedForce("COULOMB_CONSTANT * q1 * q2 * (erf(p*r)-erf(ALPHA*r)) / r")
forceGaussian.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
forceGaussian.setCutoffDistance(CUTOFF_DIST)

# forceGaussian.setForceGroup(1)
forceGaussian.addGlobalParameter("ALPHA", ALPHA)
forceGaussian.addGlobalParameter("p", p)
forceGaussian.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceGaussian.addPerParticleParameter("q")
for i in range(N_PARTICLES):
    forceGaussian.addParticle([CHARGE[i]])
systemGaussian.addForce(forceGaussian)
context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
contextGaussian = mm.Context(systemGaussian, integratorGaussian, mm.Platform.getPlatformByName('Reference'))
# context = mm.Context(system, integrator, mm.Platform.getPlatformByName('OpenCL'))
# contextGaussian = mm.Context(systemGaussian, integratorGaussian, mm.Platform.getPlatformByName('OpenCL'))
########### Calculate Directly #################
q1, q2 = CHARGE
DISTANCE_LIST, PE_PME_LIST, PE_PME_K_LIST,  PE_GAUSSIAN_LIST = [], [], [], []
context.setVelocities(np.zeros((N_PARTICLES, 3)))
contextGaussian.setVelocities(np.zeros((N_PARTICLES, 3)))
for x in arange(1.52, 4.0, 0.02):
    TwoParticles = np.asarray([[1.5, 0, 0], [x, 0, 0]])
    context.setPositions(TwoParticles)
    contextGaussian.setPositions(TwoParticles)
    state = context.getState(getEnergy=True, getPositions=True)
    state_K = context.getState(getEnergy=True, getPositions=True, groups=2)
    r = abs(state.getPositions()[0][0]-state.getPositions()[1][0])
    r1 = min(abs(state.getPositions()[0][0]-state.getPositions()[1][0]), abs(state.getPositions()[0][0]-state.getPositions()[1][0]+2.0*nanometer))
    PE_PME = state.getPotentialEnergy()
    stateGaussian = contextGaussian.getState(getEnergy=True, getPositions=True)
    PE_GAUSSIAN = stateGaussian.getPotentialEnergy()
    PE_GAUSSIAN = PE_GAUSSIAN + state_K.getPotentialEnergy()
    # PE_PME_K = state_K.getPotentialEnergy() * erf(p*r1) + erfc(ALPHA*r1) * PE_GAUSSIAN
    # PE_PME_K = state_K.getPotentialEnergy() + erfc(ALPHA*r1)*PE_GAUSSIAN - erfc(p*r1)*erf(ALPHA*r1)/erf(p*r1)/erfc(ALPHA*r1)*erfc(ALPHA*r1)*PE_GAUSSIAN
    PE_PME_K = COULOMB_CONSTANT * q1 * q2 * erf(p*r1) / r1
    DISTANCE_LIST.append(r.value_in_unit(r.unit))
    PE_GAUSSIAN_LIST.append(PE_GAUSSIAN.value_in_unit(PE_GAUSSIAN.unit))
    PE_PME_LIST.append(PE_PME.value_in_unit(PE_PME.unit))
    PE_PME_K_LIST.append(PE_PME_K.value_in_unit(PE_PME_K.unit))
    print "r: ," , r, "Gaussian ratio: ", (PE_GAUSSIAN)
###### plot ##########
plot(DISTANCE_LIST, PE_PME_LIST, 'rx-')
plot(DISTANCE_LIST, PE_GAUSSIAN_LIST, 'bx-')
plot(DISTANCE_LIST, PE_PME_K_LIST, 'yx-')
axis([0, 2.5, -300, -10])
title('PME Interaction Energy. Platform=%s' % context.getPlatform().getName())
show()





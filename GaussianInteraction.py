from math import*
import numpy as np
from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import*
from pylab import*

epsilon = 8.854187817620E-12*farad/meter

COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)

CUTOFF_DIST = 5*nanometer

#########Specify Particles###############
system = mm.System()
BOX_LENGTH = 10
system.setDefaultPeriodicBoxVectors((BOX_LENGTH, 0, 0), (0, BOX_LENGTH, 0), (0, 0, BOX_LENGTH))
N_PARTICLES = 2
MASS = 1.0
CHARGE = [1.0, -1.0]*elementary_charge
for i in range(N_PARTICLES):
    system.addParticle(MASS)
#########Specify NonBonded Force#############
force = mm.NonbondedForce()
force.setNonbondedMethod(mm.NonbondedForce.PME)
force.setReciprocalSpaceForceGroup(1)
for i in range(N_PARTICLES):
    force.addParticle(CHARGE[i], 1.0, 0.0)
system.addForce(force)
#########Specify Gaussian force##########
'''
TODO make a and b per particle parameters
TODO truncate gaussian potential with *erfc(r/ALPHA)
'''
'''
a, b = [1.0/((0.4*nanometer)**2)]*2
forceGaussian = mm.CustomNonbondedForce("COULOMB_CONSTANT * q1 * q2 * erf(p*r) / r")
# forceGaussian.setForceGroup(1)
forceGaussian.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceGaussian.addGlobalParameter("p", sqrt(a * b / (a + b)))
forceGaussian.addPerParticleParameter("q")
for i in range(N_PARTICLES):
    forceGaussian.addParticle([CHARGE[i]])
systemGaussian.addForce(forceGaussian)
'''
###########Set up integrator#################
integrator = mm.LangevinIntegrator(300, 1.0, 0.002)
# context = mm.Context(system, integrator, mm.Platform.getPlatformByName('CPU'))
context = mm.Context(system, integrator)
force.setCutoffDistance(CUTOFF_DIST)
DISTANCE_LIST, PE_NOCUTOFF_LIST, PE_DIRECT_LIST = [], [], []
for x in arange(2.0, 10.0, 0.2):
	TwoParticles = np.asarray([[1.5, 0, 0], [x, 0, 0]])
	context.setPositions(TwoParticles)
	context.setVelocities(np.zeros((N_PARTICLES, 3)))
	state = context.getState(getEnergy=True, getPositions=True)
	PE_NOCUTOFF = state.getPotentialEnergy()
	print "OpenMM PE: ", PE_NOCUTOFF
	R = abs(state.getPositions()[0][0]-state.getPositions()[1][0])
	DISTANCE_LIST.append(R.value_in_unit(R.unit))
	PE_NOCUTOFF_LIST.append(PE_NOCUTOFF.value_in_unit(PE_NOCUTOFF.unit))
plot(DISTANCE_LIST, PE_NOCUTOFF_LIST, 'bx-')
# plot(DISTANCE_LIST, PE_DIRECT_LIST)
# axis([0, 50, -40, 0])
title('Electric interaction Energy. Platform=%s' % context.getPlatform().getName())
xlabel('Displacement')
ylabel('Energy')
show()





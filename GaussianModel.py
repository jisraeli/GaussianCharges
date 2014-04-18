from math import*
import numpy as np
from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import*
from pylab import*
from sys import stdout
import simulation1

epsilon = 8.854187817620E-12*farad/meter
COULOMB_CONSTANT = (AVOGADRO_CONSTANT_NA/(4.0*pi*epsilon)).value_in_unit_system(md_unit_system)
CUTOFF_DIST = 1*nanometer
'''
Set up system and integrator
'''
system = mm.System()
N_PARTICLES = 2
MASS = 1.0 * atomic_mass_unit   
CHARGE = [1.0, -1.0] * elementary_charge
for i in range(N_PARTICLES):
    system.addParticle(MASS)
integrator = mm.CustomIntegrator(0.0001)
integrator.addPerDofVariable("x1", 0)
integrator.addUpdateContextState()
integrator.addComputePerDof("v", "v+0.5*dt*f1/m")
integrator.addComputePerDof("x", "x+dt*v")
integrator.addConstrainPositions()
integrator.addComputePerDof("v", "v+0.5*dt*f1/m+(x-x1)/dt")
integrator.addConstrainVelocities()
'''
Set up PME reciprocal and Gaussian forces in group 1 and simulate
'''
force = mm.NonbondedForce()
force.setNonbondedMethod(mm.NonbondedForce.PME)
force.setReciprocalSpaceForceGroup(1)
a, b = [1.0/((0.2*nanometer)**2)]*2
p = sqrt(a * b / (a + b))
ERROR_TOL = force.getEwaldErrorTolerance()
ALPHA = sqrt(-log(2*ERROR_TOL))/CUTOFF_DIST
forceGaussian = mm.CustomNonbondedForce("COULOMB_CONSTANT * q1 * q2 * (erf(p*r)-erf(ALPHA*r)) / r")
forceGaussian.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
forceGaussian.addGlobalParameter("p", p)
forceGaussian.addGlobalParameter("ALPHA", ALPHA)
forceGaussian.addGlobalParameter("COULOMB_CONSTANT", COULOMB_CONSTANT)
forceGaussian.addPerParticleParameter("q")
for i in range(N_PARTICLES):
    force.addParticle(CHARGE[i], 1.0, 0.0)
    forceGaussian.addParticle([CHARGE[i]])
system.addForce(force), system.addForce(forceGaussian)
TwoParticles = np.asarray([[0.5, 0, 0], [1.5, 0, 0]])
# context = mm.Context(system, integrator, mm.Platform.getPlatformByName('CPU'))
# context.setPositions(TwoParticles)
# context.setVelocities(np.zeros((N_PARTICLES, 3)))
DISTANCE_LIST, TIME_LIST ,ENERGY_LIST, POS_LIST = [], [], [], []

platform = mm.Platform.getPlatformByName('CPU')

def makeTopology(n_atoms):
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue('X', chain)
    for i in range(n_atoms):
        topology.addAtom('C', app.element.carbon, residue)
    return topology

fakeTopology = makeTopology(N_PARTICLES)
fakeTopology.setUnitCellDimensions((2,2,2))

simulation = simulation1.Simulation1(fakeTopology, system, integrator, platform)
simulation.context.setPositions(TwoParticles)
simulation.context.setVelocities(np.zeros((N_PARTICLES, 3)))
simulation.reporters.append(app.DCDReporter('trajectory.dcd', 5))
simulation.reporters.append(app.StateDataReporter(stdout, 10, step=True, 
     totalEnergy=True, separator='\t'))
simulation.step(1000)



'''
for i in range(100):
    state = context.getState(getEnergy=True, getPositions=True, groups=2)
    Positions = state.getPositions()
    POS_LIST.append(Positions)
    r = abs(state.getPositions()[0][0]-state.getPositions()[1][0])
    r1 = min(abs(state.getPositions()[0][0]-state.getPositions()[1][0]), abs(state.getPositions()[0][0]-state.getPositions()[1][0]+2.0*nanometer))
    DISTANCE_LIST.append(r.value_in_unit(r.unit))
    ENERGY = state.getPotentialEnergy()
    ENERGY_LIST.append(ENERGY.value_in_unit(ENERGY.unit))
    TIME_LIST.append(i)
    print "i, Energy", i, ENERGY
    integrator.step(1)
fakeTopology = app.Topology()
fakeTopology.setUnitCellDimensions((2,2,2))
app.pdbfile.PDBFile.writeModel(fakeTopology, POS_LIST)
'''
#plot(TIME_LIST, DISTANCE_LIST, 'rx-')
# plot(DISTANCE_LIST, PE_GAUSSIAN_LIST, 'bx-')
# plot(DISTANCE_LIST, PE_PME_K_LIST, 'yx-')
#axis([0, 50, 0, 1000000])
# title('PME Interaction Energy. Platform=%s' % context.getPlatform().getName())
#show()
'''
platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(fakeTopology, system, integrator, platform)
simulation.context.setPositions(TwoParticles)
simulation.context.setVelocities(np.zeros((N_PARTICLES, 3)))
# simulation.reporters.append(app.DCDReporter('trajectory.dcd', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 10, step=True, 
    totalEnergy=True, separator='\t'))
simulation.step(100)
'''






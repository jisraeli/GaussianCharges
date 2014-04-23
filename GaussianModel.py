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

def makeTopology(n_atoms):
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue('X', chain)
    for i in range(n_atoms):
        topology.addAtom('C', app.element.carbon, residue)
    return topology
'''
Set up -  system and integrator
'''
system = mm.System()
N_PARTICLES = 2
MASS = 1.0 * atomic_mass_unit   
CHARGE = [1.0, -1.0] * elementary_charge
for i in range(N_PARTICLES):
    system.addParticle(MASS)
integrator = mm.CustomIntegrator(0.005)
integrator.addPerDofVariable("x1", 0)
integrator.addUpdateContextState();
integrator.addComputePerDof("v", "v+0.5*dt*f1/m")
integrator.addComputePerDof("x", "x+dt*v")
integrator.addComputePerDof("x1", "x")
integrator.addConstrainPositions()
integrator.addComputePerDof("v", "v+0.5*dt*f1/m+(x-x1)/dt")
integrator.addConstrainVelocities()
'''
Set up PME reciprocal and Gaussian forces in group 1
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
platform = mm.Platform.getPlatformByName('CPU')
fakeTopology = makeTopology(N_PARTICLES)
fakeTopology.setUnitCellDimensions((2,2,2))
simulation = simulation1.Simulation1(fakeTopology, system, integrator, platform)
simulation.context.setPositions(TwoParticles)
simulation.context.setVelocities(np.zeros((N_PARTICLES, 3)))
simulation.reporters.append(app.DCDReporter('trajectory.dcd', 5))
simulation.reporters.append(app.StateDataReporter(stdout, 10, step=True, kineticEnergy=True, potentialEnergy=True,
     totalEnergy=True, separator='\t'))
simulation.step(10000)
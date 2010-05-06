import numpy 
numpy.random.seed(122222)

from amuse.support.units import nbody_system
from amuse.ext.plummer import MakePlummerModel
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from sandbox.spz.LagrangianRadii import *
from amuse.legacy.bhtree.interface import BHTreeInterface, BHTree

radius = 0.0 | nbody_system.length
N = 124
eta = 0.005
etas = 0.005
eps = 0.001 * 100.0 | nbody_system.length
dt = 0.5 | nbody_system.time
tmax = 500.0 | nbody_system.time

parts = MakePlummerModel(N).result
parts.radius = radius

#gravity=PhiGRAPE()
# BHTree code
#convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
gravity = BHTree()
gravity.parameters.epsilon_squared = eps ** 2
#gravity.parameters.timestep_parameter=eta 
#gravity.parameters.initial_timestep_parameter=etas 
#gravity.commit_parameters()
#gravity.parameters.timestep = dt / 1024.0
gravity.parameters.timestep = dt / 64.0
gravity.parameters.opening_angle = 0.5 | units.none
gravity.commit_parameters()
#gravity.commit_particles()
gravity.particles.add_particles(parts)


time = 0 | nbody_system.time
#e0=gravity.particles.kinetic_energy() + \
 #gravity.particles.potential_energy(smoothing_length_squared=eps**2,G=nbody_system.G)
#e0=0.0#gravity.kinetic_energy + gravity.potential_energy

# array for storing results
MassFraction = [0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] \
    | units.none
res = []
resCounter = 0

# writing to the file
DAT = open('core_collapse_dat', 'w')

e0 = None
while time < tmax:
  
  time = time + dt
  gravity.evolve_model(time)
  
  ek = gravity.kinetic_energy
  ep = gravity.potential_energy
  if (e0 is not None):
	print "time,T/V,err_E:", time, (ek / ep).number, ((ep + ek - e0) / e0).number
  else:
	e0 = gravity.kinetic_energy + gravity.potential_energy
  
  # calculate lagrangian radii
  Lagrad = LagrangianRadii(parts, verbose=1)
  for i in Lagrad:
	res.append([])
	res[resCounter].append(time.number)
	res[resCounter].append(MassFraction[resCounter].number)
	res[resCounter].append(i.number)
	resCounter = resCounter + 1
	DAT.write(str(res[resCounter - 1][0]) + " " + str(res[resCounter - 1][1]) + " " + str(res[resCounter - 1][2]) + " \n")
  #LagrangianRadii(gravity.particles, verbose=1)
  
  # search for bound binaries
  for a in gravity.particles:
	  for b in gravity.particles:
		  # check whether star a and b are binary
		  Ek = 0.5 * (a.)
  
  break # TODO DEBUG
  
# print results
for a in res:
	print a

DAT.close()
del gravity

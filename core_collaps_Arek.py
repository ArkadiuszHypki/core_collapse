import numpy 
numpy.random.seed(122222)

from amuse.support.units import nbody_system
from amuse.ext.plummer import MakePlummerModel
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from sandbox.spz.LagrangianRadii import *
from amuse.legacy.bhtree.interface import BHTreeInterface, BHTree

DEBUG = 1
radius = 0.0 | nbody_system.length
N = 124
eta = 0.005
etas = 0.005
# smoothing value
#eps = 0.001 * 100.0 | nbody_system.length
eps = 0.001 | nbody_system.length
dt = 0.5 | nbody_system.time
tmax = 500.0 | nbody_system.time
isBinary = numpy.zeros(N)
epsTolerance = 1000 * eps.number
  
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
gravity.parameters.timestep = dt / 1024.0
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

# writing to the file
DAT = open('core_collapse_dat', 'w')

e0 = None
from_gravity_to_model = gravity.particles.new_channel_to(parts)
while time < tmax:
  res = []
  binCounts = numpy.zeros(9)
  isBinary = numpy.zeros(N)
  time = time + dt
  gravity.evolve_model(time)
  
  ek = gravity.kinetic_energy
  ep = gravity.potential_energy
  if (e0 is not None):
	print "time,T/V,err_E:", time, (ek / ep).number, ((ep + ek - e0) / e0).number
  else:
	e0 = gravity.kinetic_energy + gravity.potential_energy
	
  # save time
  res.append(time.number)
  
  # calculate lagrangian radii
  Lagrad = LagrangianRadii(parts, verbose=1)
  resCounter = 0
  for i in Lagrad:
	res.append(i[0].number)
	
	print "time",time.number,"massFrac",MassFraction[resCounter].number,"realMass",i[1].number,"radius",i[0].number
	resCounter = resCounter + 1
	
  # search for bound binaries
  i = -1
  binCount = 0
  vj = numpy.sqrt((gravity.particles.vx[0:N].number)**2 + (gravity.particles.vy[0:N].number)**2 + 
	  (gravity.particles.vz[0:N].number)**2)
  vjRadius = numpy.sqrt((gravity.particles.x[0:N].number)**2 + (gravity.particles.y[0:N].number)**2 + 
	  (gravity.particles.z[0:N].number)**2)
  pos = gravity.particles
  
  r = []
  vesci = []
  energies = []
  i=0
  for part in pos:
	  rOne = numpy.sqrt((part.x.number - pos.x[0:N].number)**2 + (part.y.number - pos.y[0:N].number)**2 + 
		  (part.z.number - pos.z[0:N].number)**2)
	  #energies = 0.5 * (gravity.particles.mass[i].number * vj[i] + gravity.particles.mass[j].number * vj[j]) -\
		#	  	(gravity.particles.mass[i].number * gravity.particles.mass[j].number / r[i][j])
	  r.append(rOne)
	  vesci.append(numpy.sqrt(2.0 * gravity.particles.mass[0:N].number/rOne))
	  
	  #vjj = vj[0:N]
	  enOne = 0.5 * (part.mass.number * vj[i] * vj[i] + pos.mass[0:N].number * vj[0:N] * vj[0:N]) -\
	  	(part.mass.number * pos.mass.number / rOne)
	  energies.append(enOne)
	  i=i+1
	  
  #print "Energies finished",energies[0]
  
  if (DEBUG == 1):
	  # sprawdzenie bo nie wierze
	  for i in range(0, 5):
		  for j in range(0, 5):
			rOne = numpy.sqrt((pos.x[i].number - pos.x[j].number)**2 + (pos.y[i].number - pos.y[j].number)**2 + 
			  (pos.z[i].number - pos.z[j].number)**2)
			enOne = 0.5 * (pos.mass[i].number * vj[i]**2 + pos.mass[j].number * vj[j]**2) -\
				(pos.mass[i].number * pos.mass[j].number / r[i][j])
			if (abs(enOne - energies[i][j]) > 1.0e-8):
				print "ERROR energy ",i,j,enOne,energies[i][j],abs(enOne - energies[i][j])
			if (abs(rOne - r[i][j]) > 1.0e-8):
				print "ERROR r",i,j,rOne,r[i][j],abs(rOne - r[i][j])
  
	  
#  for i in range(0, N-1):
#	  energies.append(numpy.zeros(N))
#	  for j in range(i + 1, N):
#		  energies[i][j] = 0.5 * (pos.mass[i].number * vj[i] + pos.mass[j].number * vj[j]) -\
#			  	(pos.mass[i].number * pos.mass[j].number / r[i][j])
#  print "Energies finished"
	  
  for i in range(0, N-1):
	  viRadius = vjRadius[i]
	  for j in range(i + 1, N):
		  if ((isBinary[i] == 0) and (isBinary[j] == 0)):
			  # check whether star a and b are binary
			  minRadius = min(viRadius, vjRadius[j])
			  if (vj[j] < vesci[i][j] and r[i][j] < epsTolerance and energies[i][j] < 0.0):
				isBinary[i] = 1
				isBinary[j] = 1
				print "binaryStar i",i,"j",j,"time",time.number,"vesci",vesci[i][j],"r",r[i][j],"vj",vj[j]
				binCount = binCount + 1
				# put binary to specific lagrangian radius
				k = 0
				for m in Lagrad:
					if (minRadius < m[0].number):
						binCounts[k] = binCounts[k] + 1
					k = k + 1
				
  
  # save binary count
  print "time",time.number,"binCount",binCount
  res.append(binCount)
  
  # save to file
  for i in xrange(len(res)):
	  DAT.write(str(res[i]) + "  ")
  for i in xrange(len(binCounts)):
	  DAT.write(str(binCounts[i]) + "  ")
  DAT.write("\n")
    
# print results
print "res",res

DAT.close()
del gravity

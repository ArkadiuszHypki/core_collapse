import numpy
numpy.random.seed(122222)
import operator
from amuse.support.units import nbody_system
from amuse.ext.plummer import MakePlummerModel
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.support.units import units
from amuse.support.io import store
from amuse.support.io import write_set_to_file
import pickle
import os

radius = 0.0 |  nbody_system.length
N=1024
eta=0.005
etas=0.005
eps=0.001 | nbody_system.length
dt=0.5 | nbody_system.time
tmax=500 | nbody_system.time


parts=MakePlummerModel(N).result
#Mass segregation
N2 = 82
x = 1.25
N1 = N-N2
m1 = 1/(N1+x*N2)
m2 = x*m1
M1 = numpy.zeros((N1,1)) + (m1)
M2 = numpy.zeros((N2,1)) + (m2)
M = numpy.concatenate((M1,M2))
parts.mass = nbody_system.mass.new_quantity(numpy.hstack(M))
parts.radius=radius

gravity=PhiGRAPE()#(mode="gpu")
gravity.parameters.epsilon_squared = eps**2
gravity.particles.add_particles(parts)
gravity.parameters.timestep_parameter=eta 
gravity.parameters.initial_timestep_parameter=etas 
time=0 | nbody_system.time
e0=gravity.particles.kinetic_energy() + \
 gravity.particles.potential_energy(smoothing_length_squared=eps**2,G=nbody_system.G)

MassFraction = [0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] \
    | units.none
def distance_sq(stars) :
    return stars.x**2  + stars.y**2  +stars.z**2


def LagrangianRadii(stars, verbose = 0, massf = MassFraction) :

    com = stars.center_of_mass()
    stars.position = stars.position - com

    vcom = stars.center_of_mass_velocity()
    stars.velocity = stars.velocity - vcom

    n_stars = len(stars)
    # Next line is a potential performance bottleneck, becasue
    # the for loop is not in numpy but explicit
    # old but slow: d2 = numpy.array([distance_sq(star) for star in stars])
    d2 = distance_sq(stars)
    m = stars.mass/stars.mass.sum()
    d2m = zip(d2, m)
    d2m.sort(key = operator.itemgetter(0))

    iL = 0
    mt = 0 | units.none
    Lagrad = []
    for d2i,mi in d2m :
        mt += mi
        while mt >= massf[iL] :
            Lagrad.append(d2i.sqrt())
            # Previous line is preferable above here below, beacuse
            #  the former preserves the units
            #  Lagrad.append(sqrt(d2[ni].value_in(nbody_system.length**2)))
            if verbose==1 :
                print "Lagrangian Radius M= ", mt, \
                      "(iL=", iL, ") at d= ", Lagrad[-1]
            iL += 1
            if iL >= len(massf) :
                break
    return Lagrad

output_file = "nbody.hdf5"
if os.path.exists(output_file):
   os.remove(output_file)
storage = store.StoreHDF(output_file)

from_gravity_to_model = gravity.particles.new_channel_to(parts)
#L = []
#output = open('data.pkl', 'wb')
while time < tmax:
  ek=gravity.kinetic_energy
  ep=gravity.potential_energy
  print "time,T/V,err_E:", time,(ek/ep).number, ((ep+ek-e0)/e0).number
#  L.append(LagrangianRadii(parts)[-1])
  from_gravity_to_model.copy()    
  parts.savepoint(time)
  storage.store(parts.previous_state())
  time=time+dt
  gravity.evolve_model(time)
#  write_set_to_file(parts,'cc_nb/cc_nb_'+str(time.number/dt.number)+'.txt')
  #parts_temp = gravity.particles.copy()
  #print parts_temp
#  pickle.dump(adaas,output)

storage.close()
#output.close()

del gravity

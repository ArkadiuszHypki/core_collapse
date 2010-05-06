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
import sys
from amuse.legacy.fi.interface import Fi
from amuse.support.data import core
from amuse.legacy.support import channel
def main():
    N=10
    #convert_nbody = nbody_system.nbody_to_si(1.0*N | units.MSun, 149.5e3 | units.AU)
    
    radius = 0.0 |  nbody_system.length
    #radius = convert_nbody.to_si(radius)
    eta=0.05 | nbody_system.time
    #eta = convert_nbody.to_si(eta) 
    etas=0.005
    eps=0.19 | nbody_system.length
    #eps = convert_nbody.to_si(eps) 
    dt=0.5 | nbody_system.time
    #dt = convert_nbody.to_si(dt)
    tmax=10 | nbody_system.time
    #tmax = convert_nbody.to_si(tmax)
    
    parts=MakePlummerModel(N).result
    parts.radius=radius
    #gravity = Fi(convert_nbody)
    gravity = Fi()
    gravity.initialize_code()
    gravity.parameters.epsilon_squared = eps**2
    gravity.parameters.timestep = eta
    #gravity.parameters.initial_timestep_parameter=etas 
    gravity.commit_parameters()
    gravity.particles.add_particles(parts)
    time=0 | nbody_system.time
    #time = convert_nbody.to_si(time)
    e0=gravity.particles.kinetic_energy() + \
     gravity.particles.potential_energy()
    #e0=gravity.particles.kinetic_energy() + \
    # gravity.particles.potential_energy(smoothing_length_squared=eps**2,G=nbody_system.G)
    
    output_file = "Fi.hdf5"
    if os.path.exists(output_file):
       os.remove(output_file)
    storage = store.StoreHDF(output_file)
    
    #from_gravity_to_model = gravity.particles.new_channel_to(parts)
    #L = []
    #output = open('data.pkl', 'wb')
    parts = gravity.particles.copy()
    while time < tmax:
      ek = gravity.kinetic_energy
      ep = gravity.potential_energy
      print "time,T/V,err_E:", time,(ek/ep).number, ((ep+ek-e0)/e0).number
      #print "time,T/V,err_E:", time, ep
    #  L.append(LagrangianRadii(parts)[-1])
      #gravity.particles.new_channel_to(parts)
      #parts = gravity.particles.copy()
      parts.savepoint(time)
      storage.store(parts.previous_state())
      time=time+dt
      gravity.evolve_model(time)
    #  write_set_to_file(parts,'cc_nb/cc_nb_'+str(time.number/dt.number)+'.txt')
      #parts_temp = gravity.particles.copy()
      #print parts_temp
    #  pickle.dump(adaas,output)
    parts.savepoint(time)
    storage.store(parts.previous_state())
    storage.close()
    #output.close()
    
    del gravity
    del parts
if __name__ == '__main__':
	main()
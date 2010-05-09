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
    N=1024
    N2 = 0
    x = 1.25

    #convert_nbody = nbody_system.nbody_to_si(1.0*N | units.MSun, 149.5e3 | units.AU)
    
    radius = 0.0 |  nbody_system.length
    #radius = convert_nbody.to_si(radius)
    eta=1/32. | nbody_system.time
    #eta = convert_nbody.to_si(eta) 
    eps=0.001 | nbody_system.length
    #eps = convert_nbody.to_si(eps) 
    dt=0.5 | nbody_system.time
    #dt = convert_nbody.to_si(dt)
    tmax=350 | nbody_system.time
    #tmax = convert_nbody.to_si(tmax)
    
    parts=MakePlummerModel(N).result
    parts.radius=radius
    #Mass segregation
    if N2 >0:
     N1 = N-N2
     m1 = 1/(N1+x*N2)
     m2 = x*m1
     M1 = numpy.zeros((N1,1)) + (m1)
     M2 = numpy.zeros((N2,1)) + (m2)
     M = numpy.concatenate((M1,M2))
     parts.mass = nbody_system.mass.new_quantity(numpy.hstack(M))
    #gravity = Fi(convert_nbody)
    gravity = Fi()
    gravity.initialize_code()
    gravity.parameters.epsilon_squared = eps**2
    gravity.parameters.timestep = eta
    gravity.legacy_interface.set_gdgtol(0.001)
    gravity.legacy_interface.set_gdgop(0)
    gravity.legacy_interface.set_usequad(0)
    gravity.legacy_interface.set_acc_tstp(0)
    gravity.legacy_interface.set_bh_tol(0.5)
    gravity.legacy_interface.set_tstpcr2(0.025)
    gravity.legacy_interface.set_sqrttstp(0)
    gravity.legacy_interface.set_tstepcrit(0.25)
    #gravity.parameters.initial_timestep_parameter=etas 
    gravity.commit_parameters()
    gravity.particles.add_particles(parts)
    time=0 | nbody_system.time
    #time = convert_nbody.to_si(time)
    gravity.synchronize_model()
    #e0=gravity.particles.kinetic_energy() + \
    # gravity.particles.potential_energy()
#    e0=gravity.particles.kinetic_energy() + \     gravity.particles.potential_energy(smoothing_length_squared=eps**2,G=nbody_system.G)
    e0=gravity.kinetic_energy + \
     gravity.potential_energy
    
    output_file = "Fi.hdf5"
    if os.path.exists(output_file):
       os.remove(output_file)
    storage = store.StoreHDF(output_file)
    
    #from_gravity_to_model = gravity.particles.new_channel_to(parts)
    #L = []
    #output = open('data.pkl', 'wb')
    parts = gravity.particles.copy()
#    ek = gravity.kinetic_energy
#    ep = gravity.potential_energy
#    e0 = ek + ep
    while time < tmax:
      ek = gravity.kinetic_energy
      ep = gravity.potential_energy
      print "time,T/V,err_E:", time,(ek/ep).number, ((ep+ek-e0)/e0).number
      #print "time,T/V,err_E:", time, ep
    #  L.append(LagrangianRadii(parts)[-1])
      gravity.particles.new_channel_to(parts).copy()
      #parts = gravity.particles.copy()
      #parts = gravity.particles.copy()
      parts.savepoint(time)
      storage.store(parts.previous_state())
      time=time+dt
      gravity.evolve_model(time)
      gravity.synchronize_model()
    #  write_set_to_file(parts,'cc_nb/cc_nb_'+str(time.number/dt.number)+'.txt')
      #parts_temp = gravity.particles.copy()
      #print parts_temp
    #  pickle.dump(adaas,output)
#    parts.savepoint(time)
#    storage.store(parts.previous_state())
    storage.close()
    #output.close()
    
    del gravity
    del parts
if __name__ == '__main__':
	main()

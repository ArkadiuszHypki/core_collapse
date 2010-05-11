from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.io import store
from LagrangianRadii import LagrangianRadii
from matplotlib import pyplot
import os
from  amuse.support.io import read_set_from_file
import numpy
def get_Stars(filename):
        storage = store.StoreHDF(filename)
        stars = storage.load()
        #storage.close()
        Stars = []
        while stars:
            Stars.append(stars)
            stars = stars.previous_state()
        return Stars[::-1]

def get_LagrangianRadii(Stars):
    L = []
    for stars in Stars:
        L.append(LagrangianRadii(stars))
    return numpy.array(L)
    
def get_time(Stars):
    t = []
    for stars in Stars:
        t.append(stars.get_timestamp())
    return numpy.array(t.number)
 
def load_dir(dir,format='tsf'):
    print "H" 
    files = numpy.array(os.listdir(dir))
    t = numpy.array([f[10:-4] for f in files],'f')
    #sort
    files = files[numpy.argsort(t)]
    print files
    Stars = []
    for f in files:
     if (f[-3:]==format):
         Stars.append(read_set_from_file('./octgrav/'+f,format=format))
    return Stars

def Get_Mass_seg(Stars):
    N2 = sum(Stars[0].mass.number == Stars[0].mass.number.max())
    StarsM1 = [s[:-N2] for s in Stars]
    StarsM2 = [s[-N2:] for s in Stars]
    return StarsM1, StarsM2

def remove(Stars,pr):
    S = []
    for s in Stars:
        r = numpy.sum((s.position.number - s.center_of_mass().number)**2,1)
        index = numpy.argsort(r)
        npr = int(len(Stars[0])*pr)
        keys = s._get_keys()[index[-npr:]]
        st = s.copy()
        st._remove_particles(keys)
        S.append(st)
    return S

def get_mid_r(Stars):
   r = []
   for s in Stars:
        rt=numpy.sum((s.position.number - s.center_of_mass().number)**2,1)
        r.append(numpy.sort(rt)[int(rt.size/2.)])
   return numpy.sqrt(r)
   
def get_mean_r(Stars):
   r = []
   for s in Stars:
        rt=numpy.sum((s.position.number - s.center_of_mass().number)**2,1)
        r.append(rt.mean())
   return numpy.sqrt(r)

def plot_stars(S1,S2):
    pyplot.hold(False)
    for s1,s2 in zip(S1,S2):
        pyplot.plot(s1.x.number,s1.y.number,'x',s2.x.number,s2.y.number,'o')
        
def plot_L(L,t):
    pyplot.semilogy(t/2**(1.5),L[:,[4,5,6]])
    pyplot.xlabel('$t/t_D')
    pyplot.ylabel('L [0.1,0.25.0.5]')
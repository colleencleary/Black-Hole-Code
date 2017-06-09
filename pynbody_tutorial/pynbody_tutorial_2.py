# coding: utf-8
import pynbody
import numpy as np
f=pynbody.load("testdata/test_g2_snap")
len(f)
f.families()
len(f.star)
f.properties
f.properties['a']
f.keys()
f.loadable_keys()
f['pos']
f.gas.loadable_keys()
f.gas['rho']
f['double_mass']=2*f['mass']
@pynbody.derived_array
def thricethemass(sim) :
    return sim['mass']*3
f['thricethemass']
f['mass'][0]=1
f['thricethemass']
f['mass'].units
f['pos'].units
f['pos'].in_units('Mpc')
f['pos'].units
f['pos']
f['pos'].in_units('Mpc')
f.physical_units()
f['pos']
f['vel']
f.properties
f.loadable_keys()
array=pynbody.array.SimArray(np.random.rand(10))
array.sim=f
array.units='Mpc a'
array.in_units('kpc')
every_tenth=f[::10]
len(every_tenth)
every_tenth['pos'][1]
every_tenth['pos'][1]=[1,2,3]
every_tenth['pos'][1]
f['pos'][10]
f_slab=f[(f['x']>1000)&(f['x']<2000)]
f_slab['x'].min()
f_slab['x'].max()
f['x'].max()
f['x'].min()
from pynbody.filt import *
f_sphere=f[Sphere('10 kpc')]

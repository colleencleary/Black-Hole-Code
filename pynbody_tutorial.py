# coding: utf-8
import pynbody
s = pynbody.load('testdata/g15784.lr.01024.gz')
h=s.halos()
print h1
h1=h[1]
print('gas=%e)
print('gas=%e, dark matter=%e, stars=%e'%(len(h1.gas), len(h1.dark), len(h1.star)))
s
pynbody.analysis.halo.center(h1,mode='hyb')
my_h5_transform=pynbody.analysis.halo.center(h[5],mode='hyb', move_all=False)
my_h5_transform.revert()
with pynbody.analysis.halo.center(h[5], mode='hyb'): print(h[5]['pos'][0])
print(h[5]['pos'][0])
s=pynbody.load('testdata/g15784.lr.01024.gz'); h1=s.halos()[1];
cen_hyb=pynbody.analysis.halo.center(h1, mode='hyb', retcen=True)
cen_pot=pynbody.analysis.halo.center(h1, mode='pot', retcen=True)
print(cen_hyb)
print(cen_pot)
s['pos']-=cen_hyb
s.physical units()
s.physical_units()
pynbody.plot.image(h1.g, width=100, cmap='Blues');
pynbody.plot.image(s.d[pynbody.filt.Sphere('10 Mpc')], width='10 Mpc', units='Msol kpc^-2',cmap='Greys');

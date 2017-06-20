import pynbody
from pylab import *
import matplotlib.pylab as plt

s = pynbody.load('testdata/g15784.lr.01024.gz')
h=s.halos()
h1=h[1]
#print('gas=%e, dark matter=%e, stars=%e'%(len(h1.gas), len(h1.dark), len(h1.star)))
pynbody.analysis.halo.center(h1,mode='hyb')
my_h5_transform=pynbody.analysis.halo.center(h[5],mode='hyb', move_all=False)
my_h5_transform.revert()
#with pynbody.analysis.halo.center(h[5], mode='hyb'): print(h[5]['pos'][0])
#print(h[5]['pos'][0])
s=pynbody.load('testdata/g15784.lr.01024.gz'); h1=s.halos()[1];
cen_hyb=pynbody.analysis.halo.center(h1, mode='hyb', retcen=True)
cen_pot=pynbody.analysis.halo.center(h1, mode='pot', retcen=True)
#print(cen_hyb)
#print(cen_pot)
s['pos']-=cen_hyb
s.physical_units()
pynbody.plot.image(h1.g, width=100, cmap='Blues');
pynbody.plot.image(s.d[pynbody.filt.Sphere('10 Mpc')], width='10 Mpc', units='Msol kpc^-2',cmap='Greys');
pynbody.analysis.angmom.sideon(h1, cen=(0,0,0))
pynbody.plot.image(h1.g, width=100, cmap='Blues');
s.rotate_x(90)
ps=pynbody.analysis.profile.Profile(h1.s,min=0.1, max=50, type='log')
plt.clf()
plt.plot(ps['rbins'],ps['density']);
plt.semilogy();
plt.xlabel('$R$ [kpc]');
plt.ylabel('$\Sigma$ [M$_\odot$/kpc$^2$]');
plt.figure()
pd=pynbody.analysis.profile.Profile(h1.d, min=0.01, max=50, type='log')
pg=pynbody.analysis.profile.Profile(h1.g, min=0.01, max=50, type='log')
p=pynbody.analysis.profile.Profile(h1, min=0.01, max=50, type='log')
for prof, name in zip([p,pd,ps,pg],['total','dm','stars','gas']) : plt.plot(prof['rbins'],prof['v_circ'], label=name)
plt.xlabel('$R$ [kpc]');
plt.ylabel('$v_{circ}$ [km/s]');
plt.legend()

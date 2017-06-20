import pynbody
import matplotlib.pylab as plt
import numpy as np

def loader(Path):
    '''returns snap and halo'''
    s=pynbody.load(Path)
    s.physical_units()
    BH=s.star[pynbody.filt.BandPass('tform','-15 Gyr','0 Gyr')]    
    return (s,s.halos(),BH)

def plot_BH_pos(h, BH, grpid, width=100, attribute, cmap='Greys',save=False, snapshot=148):
    '''plots position of black holes: 
    h=snapshot (ie h148), grpid=halo/amiga id, w=width, attribute=property that the color map represents,
    units=preferred units of attribute, cmap=colormap name (string) from matplotlib, save=False (set True to save),
    snapshot=snapshot id (ie 148)'''
    
    pynbody.plot.image(h[grpid].s,qty=attribute,cmap=cmap, width=width, log=False)
    plt.plot(np.array((BH[np.where(BH['amiga.grp']==grpid)]['pos'])).T[0],np.array((BH[np.where(BH['amiga.grp']==grpid)]['pos'])).T[1],'r+')
    plt.xlim(-w/2,w/2)
    plt.ylim(-w/2,w/2)
    plt.title('h%s_h%s_%s_w%s'%(snapshot, grpid, attribute, w))
    
    if save==True:
        plt.savefig('h%s_h%s_tform_w%s.png'%(snapshot, grpid, attribute, w))
    plt.show()


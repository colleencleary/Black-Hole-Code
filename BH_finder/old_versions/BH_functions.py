import pynbody
import matplotlib.pylab as plt
import numpy as np
from sklearn.preprocessing import scale

def loader(path):
    '''returns snap and halo'''
    s=pynbody.load(path)
    s.physical_units()
    h=s.halos()
    BH=s.star[pynbody.filt.BandPass('tform','-15 Gyr','0 Gyr')]    
    return s, h, BH

def plot_BH_pos(h, BH, grpid, attribute, cmap='Greys', w=100, save=False, snapshot=148):
    '''plots position of black holes. Send halo array(h), black hole array(BH), grpid, and attribute you want to plot 
    Optional: width/preferred scope(w) (default=100), colormap(default='Greys'), save(default=False), snapshot(default=148)'''
    
    pynbody.plot.image(h[grpid].s,qty=attribute,cmap=cmap, width=w, log=False)
    plt.plot(np.array((BH[np.where(BH['amiga.grp']==grpid)]['pos'])).T[0],np.array((BH[np.where(BH['amiga.grp']==grpid)]['pos'])).T[1],'r+')
    plt.xlim(-w/2,w/2)
    plt.ylim(-w/2,w/2)
    plt.title('h%s_h%s_%s_w%s'%(snapshot, grpid, attribute, w))
    
    if save==True:
        plt.savefig('plots/h%s/h%s/h148_h%s_%s_w%s.png'%(snapshot, grpid, grpid, attribute, w))
    plt.show()

def plot_BH_pos_scaled(h, BH, grpid, snap, attribute, cmap='Greys', w=100, save=False, snapshot=148):
    '''plots position of black holes and scales attribute. Send halo array(h), black hole array(BH), grpid, and attribute you want to plot 
    Optional: width/preferred scope(w) (default=100), colormap(default='Greys'), save(default=False), snapshot(default=148)'''
    
    tform_array=scale(snap.star[attribute])
    
    pynbody.plot.image(h[grpid].s,qty=tform_array,cmap=cmap, width=w, log=False)
    plt.plot(np.array((BH[np.where(BH['amiga.grp']==grpid)]['pos'])).T[0],np.array((BH[np.where(BH['amiga.grp']==grpid)]['pos'])).T[1],'r+')
    plt.xlim(-w/2,w/2)
    plt.ylim(-w/2,w/2)
    plt.title('h%s_h%s_%s_scaled_w%s'%(snapshot, grpid, attribute, w))
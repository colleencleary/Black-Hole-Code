
# coding: utf-8

# In[1]:

import pynbody
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
get_ipython().magic(u'matplotlib inline')
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


# In[2]:

#galaxy_masses=pd.DataFrame()
#BH_halos=pd.DataFrame()
snap_id=0


# In[174]:

def loader(path):
    '''returns snapshot and halo'''
    snap=pynbody.load(path)
    snap.physical_units()
    h=snap.halos()
    return snap,h

def BH_finder(snap):
    '''returns subsnap of black holes'''
    return snap.star[pynbody.filt.BandPass('tform','-15 Gyr','0 Gyr')]

def analysis(halo, view=0):
    '''center snap on halo's center of mass (angmom gives errors if there aren't enough pixels). 
    view=1 for face on, 2 for side on, anything else will leave it as is'''
    if view==1:
        pynbody.analysis.angmom.faceon(halo)
    elif view==2:
        pynbody.analysis.angmom.sideon(halo)
    else:
        pynbody.analysis.halo.center(halo)
        
def plot_BH_pos(h, BH, attribute, cmap='Greys', w=100, lc='best',save=False, filename='plots/dg_adjTime/h%s_t%s_w%s.png'):
    '''plots position of black holes. Send halo array(h[grpid]), black hole array(BH),attribute 
    Optional: width/preferred scope(w) (default=100), colormap(default='Greys'), save(default=False)'''
    
    halo_id=int(np.unique(h['amiga.grp']))
    time=h.properties['time'].in_units('Gyr')
    
    pynbody.plot.image(h,qty=attribute,cmap=cmap, width=w, log=False)
    for i in np.where(BH['amiga.grp']==halo_id):
        plt.plot(np.array(BH[i]['pos']).T[0],np.array(BH[i]['pos']).T[1],'r+', label=BH[i]['iord'])
    plt.xlim(-w/2,w/2)
    plt.ylim(-w/2,w/2)
    plt.title('Halo %s at %s Gyr'%(halo_id, round(time, 1)))
    plt.legend(fancybox=True, loc='lc')
    plt.tight_layout()

    if save==True:
        plt.savefig(filename%(halo_id, round(time, 1), w),bbox_inches='tight', dpi=200)

        

def adjust_csv_index(df):
    df=df.rename(columns={'Unnamed: 0':'snap'})
    df=df.set_index('snap')
    
def add_mass_columns(hgrpid, grpid, snap_id):
    '''returns row with gas, stellar, and total masses for one halo tracked through time'''
    df=pd.DataFrame([[hgrpid.g['mass'].sum(), hgrpid.s['mass'].sum(), hgrpid['mass'].sum()]], index=[snap_id], columns=['H[%d] Gas'%(grpid),'H[%d] Star'%(grpid),'H[%d] Total'%(grpid)])
    df.units='Msol'
    return df
    
def mass_df(h, BHgrp, snap_id):
    df=pd.DataFrame()
    for i in BHgrp:
        df=pd.concat([df, add_mass_columns(h[i], i, snap_id)], axis=1, join='outer')
    df.units='Msol'
    return df
    
def galaxy_df(snap_id,h,BHgrp):
    '''snap_id_array[snap_id], h, np.unique(BH['amiga.grp'])'''
    galaxy_masses=pd.DataFrame()
    for i in snap_id:
        vars()['mass_%s'%(i)]=mass_df(h,BHgrp,i)
        galaxy_masses=galaxy_masses.append(vars()['mass_%s'%(i)])
    return galaxy_masses
    if save==True:
        galaxy_masses.to_csv('csvdata/%s_galaxy_masses.csv'%(snap_id))        
    
def plot_SFH(h, grpid, snap_id, c='k', save=True, mf=False, filename='plots/SFH/i_%s_SFH.png'): #%s_h[%s]_SFH.png'):
    '''plots SFH. Send halo array(h[grpid]), color(c=k) grpid, and snapshot you want to plot 
    Optional:save(default=True)
    #examples of iterative SFH:
    for i in np.unique(BH['amiga.grp']): #using black hole array
        plot_SFH(h,i,snap_id_array[snap_id], mf=False)
    for i in np.unique(iords[str(snap_id_array[snap_id])+' halo groups'].dropna()): #using iords array
        plot_SFH(h,i,snap_id_array[snap_id], mf=False)'''
    plt.clf()
    pynbody.plot.sfh(h[grpid].s,color=c, massform=mf);
    plt.title('%s_h[%s]_SFH'%(snap_id,grpid), y=1.14);
    if save==True:
        plt.savefig(filename,bbox_inches='tight', dpi=200);

def BH_halo_add(snap_id, BH, BH_halos):
    '''snap_id_array[snap_id], BH, BH_halos'''
    df=pd.DataFrame({'%s halo groups'%(snap_id):BH['amiga.grp'],'%s iords'%(snap_id):BH['iord']})
    BH_halos=pd.concat([BH_halos, df], axis=1, join='outer')
    return BH_halos

def load_snap_csvs(snap_id):
    '''snap_id_array[snap_id]. load mass and iord csvs'''
    #masses=pd.DataFrame()
    masses=pd.read_csv('csvdata/galaxy_masses.csv',index_col='snap')
    iords=pd.DataFrame()
    iords=pd.read_csv('csvdata/iords.csv',index_col='snap')
    return masses,iords
    
def get_iords_df(iords, iords_list, masses, save=False):
    for i in iords_list:
        ind=0
        vars()['iord_%s'%(i)]=pd.DataFrame(columns=['gas', 'star', 'total'])
        iords_temp=iords[i].dropna()
        vars()['iord_%s'%(i)]['halo']=iords[i]
        for j in iords_temp.index:
            if pd.notnull(iords.at[j,i]):
                temp=[masses.at[j,'H[%s] Gas'%(int(iords.at[j,i]))],masses.at[j,'H[%s] Star'%(int(iords.at[j,i]))],masses.at[j,'H[%s] Total'%(int(iords.at[j,i]))]]
            if (pd.notnull(temp)).any():
                vars()['iord_%s'%(i)].loc[int(iords_temp.index[ind]),['gas','star','total']]=temp
            else: continue
            ind=ind+1
            vars()['iord_%s'%(i)]=vars()['iord_%s'%(i)].fillna(0)
            if save==True:
                vars()['iord_%s'%(i)].to_csv('csvdata/iord_%s.csv'%(i))
                
def iord_mass_plot(iord_df, iord_id, save=False, filename='plots/galaxy_masses/iord_%s_masses.png'):
    plt.hist(range(0,len(iord_df)), weights=iord_df['total'], bins=len(iord_df),width=.85, alpha=.3, histtype='bar', label='Total', color='#666666',rwidth=1, align='mid')
    plt.hist([range(0,len(iord_df)),range(0,len(iord_df))], weights=[iord_df['gas'],iord_df['star']],bins=len(iord_df),width=.85, alpha=1, histtype='bar', stacked=True, label=['Gas','Stars'],fill=True, color=['#44AA99','#88CCEE'],rwidth=1,align='mid')
    plt.title('Masses of Halo Containing iord %s'%(iord_id));
    plt.ylabel('Mass [M$_\odot$]', fontsize='large');
    plt.legend(fancybox=True,loc='best')
    plt.yscale('log')
    if save==True:
        plt.savefig(filename%(iord_id),bbox_inches='tight')

def iord_gmass_plot(iord_df, time, bins=150, lc='best', htype= 'bar', halo_label=False, save=False, filename='plots/galaxy_masses/iord_%s_masses.png'):
    '''send individual iord dataframe, iord_id, timearray(time['time[Gyr]'])'''
    plt.clf()
    iord_id=iord_df.name
    
    #plt.hist(time['time[Gyr]'], weights=iord_df['total'], bins=bins, histtype=htype, label='Total', color='k', fill=True, linewidth=.2, alpha=.7);
    #plt.hist([time['time[Gyr]'], time['time[Gyr]']], weights=[iord_df['star'], iord_df['gas']], bins=bins, histtype=htype, label=['Stars', 'Gas'], color=['#FFEE33','#1B7837'],  stacked=True, fill=True);
    #plt.hist(time['time[Gyr]'], weights=iord_df['BH'], bins=bins, histtype=htype, label='BH', color='#662506', fill=False, linestyle='--');
    #plt.hist(time['time[Gyr]'], weights=iord_df['total'], bins=bins, histtype=htype, label='Total', color='k', fill=False, linestyle='--');
    plt.plot(time['time[Gyr]'], iord_df['total'], color='k', label='Total', linewidth=1.5)
    plt.plot(time['time[Gyr]'], iord_df['gas'], color='#1B7837', linestyle='--', label='Gas', linewidth=.8)
    plt.plot(time['time[Gyr]'], iord_df['star'], color='#FFEE33', linestyle='--', label='Stars', linewidth=.8)
    plt.plot(time['time[Gyr]'], iord_df['BH'], color='#662506', linestyle='--', label='BH', linewidth=.8)
    plt.title('Masses of Halo Containing iord %s'%(iord_id), y=1.14);
    plt.ylabel('Mass [M$_\odot$]', fontsize='large');
    plt.xlabel('Time [Gyr]', fontsize='large');
    plt.xlim(-.25,14)
    plt.legend(fancybox=True,loc=lc)
    plt.yscale('log')

    x0, x1 = plt.gca().get_xlim()
    old_axis = plt.gca()
    pz = plt.twiny()
    labelz=time_axis['z'][-6:]
    times=time_axis['time[Gyr]'][-6:]
    pz.set_xticks(times)
    pz.set_xticklabels([str(x) for x in labelz])
    pz.set_xlim(x0,x1)
    pz.set_xlabel('$z$')
    plt.sca(old_axis)

    if halo_label==True:
        prev=0
        for t in times:
            #if iord_df['halo'].str.contains('[').any():
                #iord_df['halo'] = iord_df['halo'].map(lambda x: str(x).lstrip('[').rstrip(']'))
            halo=iord_df[iord_df['age']==a]['halo'].item()
            if halo!=prev and halo!=0:
                plt.text(a, iord_df[iord_df['age']==a]['total'].values*1.1,str('%s'%(int(halo))), color='#666666')
                prev=halo
            else:continue

    plt.tight_layout()
    plt.subplots_adjust(top=1)
    if save==True:
        plt.savefig(filename%(iord_id),bbox_inches='tight', dpi=200)
        
def iord_bhacc_plot(iord_df, time, bins=300, save=False, filename='plots/accretion_rates/iord_%s_accRates.png'):
    '''send individual iord dataframe, time array'''
    plt.clf()
    iord_id=iord_df.loc[0,'iord']
    
    plt.hist(iord_df['time'], weights=iord_df['accRate'], bins=bins, histtype='step', color='k', linewidth=.9);
    plt.title('BH %s Accretion Rates'%(iord_id), y=1.14);
    plt.ylabel('Accretion Rate [M$_\odot$ yr$^{-1}$]]', fontsize='large');
    plt.xlabel('Time [Gyr]', fontsize='large');
    plt.xlim(-.25,14)
    plt.yscale('log')

    x0, x1 = plt.gca().get_xlim()
    old_axis = plt.gca()
    pz = plt.twiny()
    labelz=time_axis['z'][-6:]
    times=time_axis['time[Gyr]'][-6:]
    pz.set_xticks(times)
    pz.set_xticklabels([str(x) for x in labelz])
    pz.set_xlim(x0,x1)
    pz.set_xlabel('$z$')
    plt.sca(old_axis)

    plt.tight_layout()
    plt.subplots_adjust(top=1)
    if save==True:
        plt.savefig(filename%(iord_id),bbox_inches='tight', dpi=200)

def iord_bhmass_plot(iord_df, time, bins=300, save=False, filename='plots/distance/iord_%s_distance.png'):
    '''send individual iord dataframe, time array'''
    plt.clf()
    iord_id=iord_df.loc[0,'iord']
    
    plt.hist(iord_df['time'], weights=iord_df['mass'], bins=bins, histtype='step', color='k', linewidth=.9);
    plt.title('BH %s Mass'%(iord_id), y=1.14);
    plt.ylabel('Mass [M$_\odot$]', fontsize='large');
    plt.xlabel('Time [Gyr]', fontsize='large');
    plt.xlim(-.25,14)
    #plt.yscale('log')

    x0, x1 = plt.gca().get_xlim()
    old_axis = plt.gca()
    pz = plt.twiny()
    labelz=time_axis['z'][-6:]
    times=time_axis['time[Gyr]'][-6:]
    pz.set_xticks(times)
    pz.set_xticklabels([str(x) for x in labelz])
    pz.set_xlim(x0,x1)
    pz.set_xlabel('$z$')
    plt.sca(old_axis)

    plt.tight_layout()
    plt.subplots_adjust(top=1)
    if save==True:
        plt.savefig(filename%(iord_id), bbox_inches='tight', dpi=200)
        
def iord_bhr_plot(iord_df, time, bins=300, save=False, filename='plots/distance/iord_%s_distance.png'):
    '''send individual iord dataframe, time array'''
    plt.clf()
    iord_id=iord_df.loc[0,'iord']
    
    plt.hist(iord_df['time'], weights=iord_df['r'], bins=bins, histtype='step', color='k', linewidth=.9);
    plt.title('BH %s Radial Distance'%(iord_id), y=1.14);
    plt.ylabel('R [kpc]', fontsize='large');
    plt.xlabel('Time [Gyr]', fontsize='large');
    plt.xlim(-.25,14)
    plt.yscale('log')

    x0, x1 = plt.gca().get_xlim()
    old_axis = plt.gca()
    pz = plt.twiny()
    labelz=time_axis['z'][-6:]
    times=time_axis['time[Gyr]'][-6:]
    pz.set_xticks(times)
    pz.set_xticklabels([str(x) for x in labelz])
    pz.set_xlim(x0,x1)
    pz.set_xlabel('$z$')
    plt.sca(old_axis)

    plt.tight_layout()
    plt.subplots_adjust(top=1)
    if save==True:
        plt.savefig(filename%(iord_id),bbox_inches='tight', dpi=200)
        
def cycle_snap_id_array(snap_id_array,iords_list):
    '''loads each snap in snap_id_array and stores values to existing dataframes and saves as csvs(not as variable, so reload csv). Send snap_id_array, time dataframe, iords_list.
    Currently set to retrieve time[Gyr], z, BH[r], BH[mass]
    WARNING: This is going to be like a half hour wait or so
    Note: can send partial arrays'''
    time=pd.DataFrame(index=snap_id_array, columns=['time[Gyr]', 'z'])
    for j in snap_id_array:
        path='/data/scratch/jillian/h148/h148.cosmo50PLK.3072g3HbwK1BH.00%s/h148.cosmo50PLK.3072g3HbwK1BH.00%s'%(j,j)
        snap,h=loader(path)
        time.loc[j]['time[Gyr]']=snap.properties['time'].in_units('Gyr')
        time.loc[j]['z']=snap.properties['z']
        time.to_csv('csvdata/time.csv')
        BH=BH_finder(snap)
        for i in iords_list:
            temp=BH[np.where(BH['iord']==int(i))]['mass']
            if len(temp)>0:
                vars()['iord_%s'%(i)].loc[int(j),'BH']=temp
                vars()['iord_%s'%(i)].loc[int(j),'r[kpc]']=BH[np.where(BH['iord']==int(i))]['r']
            vars()['iord_%s'%(i)].to_csv('csvdata/iord_%s.csv'%(i), index=False)

def get_orbits(iords_list, save=False):
    '''read .orbit file and convert to dataframe. long run time!'''
    cols=['iord', 'time', 'stepNumber', 'mass','xPos', 'yPos', 'zPos','xVel','yVel','zVel','pot','accRate','delM','FB_E_released','dt','scalefactor']
    orbit_file=pd.read_csv('/data/scratch/jillian/h148/h148.cosmo50PLK.3072g3HbwK1BH.orbit', header=None, sep=' ', index_col=False,names=cols)
    matches=np.isin(np.array(orbit_file['iord']), np.array(iords_list))
    bhorbit=orbit_file[matches]
    if save==True:
        bhorbit.to_csv('csvdata/bhorbit.csv', index=False)

def get_iorbits(iord_id, bhorbit, save=False, filename='csvdata/i%s_orbit.csv'):
    '''return orbit dataframe for specified iord. long run time!'''
    iord_id=int(iord_id)
    mask=np.isin(np.array(bhorbit['iord']), iord_id)
    i_orbit=bhorbit[mask]
    #mask2=np.isin(temp['stepNumber'], map(int, snap_id_array))
    i_orbit.loc[:,'accRate']=(i_orbit['accRate']*1.84793e16)/38759428183.8 #M_sol/yr
    i_orbit.loc[:,'time']=(i_orbit['time']*38759428183.8)/10**9 #Gyr
    i_orbit.loc[:,'mass']=i_orbit['mass']*1.84793e16
    i_orbit.loc[:,'xPos']=i_orbit['xPos']*50000 #kpc
    i_orbit.loc[:,'yPos']=i_orbit['yPos']*50000 #kpc
    i_orbit.loc[:,'zPos']=i_orbit['zPos']*50000 #kpc
    if save==True:
        i_orbit.to_csv(filename%(iord_id), index=False)
    return i_orbit

def getr(i_orbit):
    #vars()['i%s_orbit'%(i)]=vars()['i%s_orbit'%(i)].drop('Unnamed: 0.1', axis=1)
    combined=np.vstack((i_orbit['xPos'], i_orbit['yPos'], i_orbit['zPos'])).T
    r=np.sqrt((combined ** 2).sum(axis=1))
    i_orbit['r']=pd.Series(r)
    return i_orbit

def plot_all(h, grpid, time, iord_data, iord_df, bins=200, save=False, filename='plots/dg_adjTime/i_%s_all.png'):
    #from sklearn import preprocessing
    from sklearn.preprocessing import normalize
    
    #h[grpid].s['tform']=normalize(np.arrayh[grpid].s['tform'].reshape(-1, 1))
    totalm_norm=normalize(iord_data['total'].values.reshape(-1, 1))
    gasm_norm=normalize(iord_data['gas'].values.reshape(-1, 1))
    starm_norm=normalize(iord_data['star'].values.reshape(-1, 1))
    bhsnapm_norm=normalize(iord_data['BH'].values.reshape(-1, 1))
    #accR_norm=normalize(iord_df['accRate'].values.reshape(-1, 1))
    #bhm_norm=normalize(iord_df['mass'].values.reshape(-1, 1))
    #bhr_norm=normalize(iord_df['r'].values.reshape(-1, 1))
    
    plt.clf()
    iord_id=iord_df.loc[0,'iord']
        
    #pynbody.plot.sfh(h[grpid].s, massform=False, label='SFR')
    plt.plot(time['time[Gyr]'], totalm_norm, label='Total Mass', linewidth=.9)
    plt.plot(time['time[Gyr]'], gasm_norm, linestyle='--', label='Gas Mass', linewidth=.8)
    plt.plot(time['time[Gyr]'], starm_norm, linestyle='--', label='Star Mass', linewidth=.9)
    plt.plot(time['time[Gyr]'], bhsnapm_norm, color='#662506', linestyle='--', label='BH Mass', linewidth=.8)
    plt.hist(iord_df['time'], weights=iord_df['accRate'], bins=bins, histtype='step', label='BH Accretion Rate', linewidth=.9, normed=1)
    plt.hist(iord_df['time'], weights=iord_df['mass'], bins=bins, histtype='step', label='BH Mass', linewidth=.9, normed=1);
    plt.hist(iord_df['time'], weights= iord_df['r'], bins=bins, histtype='step', label='BH Radial Distance', linewidth=.9, normed=1);
    
    #
    #pynbody.plot.sfh(h[grpid].s, label='SFR') 
    #is made of:
    #trange=[h[grpid].s['tform'].in_units("Gyr").min(), h[grpid].s['tform'].in_units("Gyr").max()]
    #trangefilt=filt.And(filt.HighPass('tform', str(trange[0]) + ' Gyr'),filt.LowPass('tform', str(trange[1]) + ' Gyr'))
    #tforms=h[grpid].s[trangefilt]['tform'].in_units('Gyr')
    #weight = simstars[trangefilt]['mass'].in_units('Msol') * binnorm
    #plt.hist(tforms, weights=weight, bins=bins, histtype='step', label='SFR')
    
    plt.xlabel('Time [Gyr]', fontsize='large');
   # plt.xlim(-.25,14)
    plt.yscale('log')

    x0, x1 = plt.gca().get_xlim()
    old_axis = plt.gca()
    pz = plt.twiny()
    labelz=time_axis['z'][-6:]
    times=time_axis['time[Gyr]'][-6:]
    pz.set_xticks(times)
    pz.set_xticklabels([str(x) for x in labelz])
    pz.set_xlim(x0,x1)
    pz.set_xlabel('$z$')
    plt.sca(old_axis)

    plt.legend(fancybox=True,)
    plt.tight_layout()
    plt.subplots_adjust(top=1)
    if save==True:
        plt.savefig(filename%(iord_id),bbox_inches='tight', dpi=200)
        
#handy:
#for i in iords_list:
#    vars()['i_%s_data'%(i)].name=i
#    iord_gmass_plot(vars()['i_%s_data'%(i)],time, lc=7, save=True)
#    iord_bhacc_plot(vars()['i%s_orbit'%(i)], time, bins=200,save=True)
#    iord_bhmass_plot(vars()['i%s_orbit'%(i)], time, bins=200,save=True)
#    iord_bhr_plot(vars()['i%s_orbit'%(i)], time, bins=200,save=True)

#maybe later?
#snap_id=snap_id-1
#snap_id_array[snap_id]
#path_higherz='/data/scratch/jillian/h148/h148.cosmo50PLK.3072g3HbwK1BH.00%s/h148.cosmo50PLK.3072g3HbwK1BH.00%s'%(snap_id_array[snap_id],snap_id_array[snap_id])
#snap_higherz,h_higherz=loader(path_higherz)
#mask=np.isin(np.array(snap['iord']), iords_list)
#BH_snap=pynbody.snapshot.IndexedSubSnap(snap, np.where(mask))
#b=pynbody.bridge.OrderBridge(BH_snap, snap_higherz)
#cat=b.fuzzy_match_catalog()

#to adjust time:
#for i in iords_list:
#    vars()['adj_i%s_orbit'%(i)]=pd.read_csv('csvdata/i%s_orbit.csv'%(i))
#    vars()['adj_i%s_orbit'%(i)].loc[:,'time']=vars()['adj_i%s_orbit'%(i)]['time']*1.0774702
#    vars()['adj_i%s_orbit'%(i)].to_csv('csvdata/adj_i%s_orbit.csv'%(i), index=False)


# In[4]:

#array of last four digits for each snap
snap_id_array=['0139','0225','0275','0640','0974','1024','1269','1280','1408','1740','2048','2088','2432','2688','2816','2944','3072','3195','3200','3328','3456','3584','3606','3712','3840','3968','4096']


# In[5]:

#snap id index (negative values start from end)
snap_id=snap_id-1
snap_id_array[snap_id]


# In[6]:

masses,iords=load_snap_csvs(snap_id_array[snap_id])
iords_list=iords.columns.tolist()
iords_list=map(int, iords_list)
time=pd.read_csv('csvdata/time.csv', index_col=0)
global time_axis
time_axis=pd.read_csv('csvdata/time_axis.csv', index_col=0)


# In[7]:

for i in iords_list:
    vars()['i_%s_data'%(i)]=pd.read_csv('csvdata/iord_%s.csv'%(i),index_col='snap')
    vars()['i%s_orbit'%(i)]=pd.read_csv('csvdata/i%s_orbit.csv'%(i))
    vars()['adj_i%s_orbit'%(i)]=pd.read_csv('csvdata/adj_i%s_orbit.csv'%(i))


# In[217]:

#set path
path='/data/scratch/jillian/h148/h148.cosmo50PLK.3072g3HbwK1BH.00%s/h148.cosmo50PLK.3072g3HbwK1BH.00%s'%(snap_id_array[snap_id],snap_id_array[snap_id])
#test for working at home:
#path='/Users/Owner/Black_Hole_Research/old/pynbody_tutorial/testdata/g15784.lr.01024.gz'
#returns snapshot and halos in physical units (takes a couple of minutes)
snap,h=loader(path)


# In[221]:

BH=BH_finder(snap)
BH_arm=BH[pynbody.filt.BandPass('amiga.grp',2,15)]


# In[105]:

for i in BH_arm['iord']:
    vars()['i_%s_data'%(i)].name=i
    for j in BH_arm['amiga.grp']:
        plot_SFH(h,j,snap_id_array[snap_id], mf=False, save=True, filename='plots/dg_adjTime/i_%s_SFH.png'%(i))
    iord_gmass_plot(vars()['i_%s_data'%(i)],time, lc=7, save=True, filename='plots/dg_adjTime/i_%s_gmass.png')
    iord_bhacc_plot(vars()['adj_i%s_orbit'%(i)], time, bins=200,save=True, filename='plots/dg_adjTime/i_%s_bhacc.png')
    iord_bhmass_plot(vars()['adj_i%s_orbit'%(i)], time, bins=200,save=True, filename='plots/dg_adjTime/i_%s_bhmass.png')
    iord_bhr_plot(vars()['adj_i%s_orbit'%(i)], time, bins=200,save=True, filename='plots/dg_adjTime/i_%s_bhr.png')


# In[237]:

path_higherz='/data/scratch/jillian/h148/h148.cosmo50PLK.3072g3HbwK1BH.00%s/h148.cosmo50PLK.3072g3HbwK1BH.00%s'%(2944,2944)
snap_higherz,h_higherz=loader(path_higherz)
#BH=BH_finder(snap)


# In[238]:

halo_id=BH[np.where(BH['iord']==101863883)]['amiga.grp']


# In[239]:

analysis(h[halo_id],view=1)


# In[251]:

plot_BH_pos(h[halo_id].s, BH,'tform', w=20, save=True)


# In[247]:

h1_at_low_z = b(h_higherz[4])


# In[224]:

b=pynbody.bridge.OrderBridge(snap_3072,snap_higherz)


# In[225]:

cat=b.fuzzy_match_catalog()


# In[232]:

snap_3072.halos()[1]


# In[235]:

cat


# In[222]:

for j in BH_arm['amiga.grp']:
        plot_SFH(h,j,snap_id_array[snap_id], mf=False, save=True, filename='plots/dg_adjTime/i_%s_SFH.png'%(i))
        plot_SFH(h, grpid, snap_id,


# In[223]:

for i in np.unique(BH_arm['amiga.grp']): #using black hole array
        plot_SFH(h,i,snap_id_array[snap_id], mf=False, save=True, filename='plots/dg_adjTime/i_%s_SFH.png'%(i))


# In[ ]:




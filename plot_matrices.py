#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 17:53:06 2018

plots the matrices from Fig. 6C
and saves the data to mat files

@author: kwilmes
"""

from pypet import Trajectory
import matplotlib.pyplot as plt
import colormaps as cmaps
import seaborn as sns; sns.set(color_codes=True)
import pandas as pd
import numpy as np
import matplotlib
import scipy.io

plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.set_cmap(cmaps.viridis)

def plot_NC(rsc_frame, varied_param, other_param, identifier):
    plt.figure()
    # set width of bar
    barWidth = 0.25
     
    # set height of bar
    bars1 = rsc_frame[(rsc_frame[other_param] == 1)&(rsc_frame[varied_param] == 0)].values[0][2:]
    bars2 = rsc_frame[(rsc_frame[other_param] == 1)&(rsc_frame[varied_param] == 1)].values[0][2:]
    bars3 = rsc_frame[(rsc_frame[other_param] == 1)&(rsc_frame[varied_param] == 2)].values[0][2:]
     
    print('rsc_frame')
    #print(rsc_frame[(rsc_frame[other_param] == 1)&(rsc_frame[varied_param] == 2)])
    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
     
    # Make the plot
    plt.bar(r1, bars1, color=cmaps.viridis(.8), width=barWidth, edgecolor='white', label='0')
    plt.bar(r2, bars2, color=cmaps.viridis(.5), width=barWidth, edgecolor='white', label='1')
    plt.bar(r3, bars3, color=cmaps.viridis(.2), width=barWidth, edgecolor='white', label='2')
     
    # Add xticks on the middle of the group bars
    plt.xlabel('Pair', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bars1))], ['PYRPYR','SSTSST','PVPV','VIPVIP','PYRSST','PYRPV','PYRVIP','SSTPV','SSTVIP','VIPPV'])
    plt.ylabel('r_sc', fontweight='bold')
     
    # Create legend & Show graphic
    plt.legend(title=varied_param)
    plt.savefig('%s/NC_change_with_%s_%s.eps'%(savepath,varied_param,identifier))#, bbox_extra_artists=(lgd,), bbox_inches='tight',format='pdf', transparent=True) 


def plot_var(var,varname,varied_param, other_param, identifier, var2=None, var2name =''):
    plt.figure(figsize=(2.5,2.5))
    ax = plt.subplot(2,1,1)
    ##Let's iterate through the columns and plot the different firing rates :
    #for param_value in SI[other_param]:
    print(var[varied_param])
    print(var[varname])
    ax = var.pivot(index=other_param,columns=varied_param,values=varname).plot(colormap=cmaps.viridis)#color=cmaps.viridis(param_value/2.0))
    if var2 is not None:
        var2.pivot(index=other_param,columns=varied_param,values=var2name).plot(colormap=cmaps.magma,ax=ax)#color=cmaps.viridis(param_value/2.0))
    #ax.set_yscale('log')
    #lgd = plt.legend(title = other_param, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False, framealpha=1)
    #plt.ylim(0,4.5)
    plt.tight_layout()
    plt.xlabel(other_param)#, fontsize=20)
    plt.ylabel(varname)#, fontsize=20)
    plt.tight_layout
    plt.savefig('%s/%s%s_with_%s_and_%s_%s.pdf'%(savepath,varname,var2name,varied_param, other_param, identifier))#, bbox_extra_artists=(lgd,), bbox_inches='tight',format='pdf', transparent=True) 

def tsplot(ax,mean,std,**kw):
    x = np.arange(mean.shape[0])
    cis = (mean - std, mean + std)
    #ax.grid(True)#, which='major', color='w', linewidth=1.0)
    ax.fill_between(x,cis[0],cis[1],alpha=0.2,**kw)
    ax.plot(x,mean,**kw, lw=2)
    ax.margins(x=0)


from mpl_toolkits.axes_grid1 import AxesGrid

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


plot = False


for identifier in ['multiplicative','additive']:

    if identifier == 'multiplicative':
        modeltype = 'mult'
    elif identifier == 'additive':
        modeltype = 'add'
    else:
        modeltype = 'mult'

    savepath = './'
    
    # This time we don't need an environment since we just going to look
    # at data in the trajectory
    traj = Trajectory('%s'%identifier, add_time=False)
    
    # Let's load the trajectory from the file
    # Only load the parameters, we will load the results on the fly as we need them
    traj.f_load(filename='%s%s.hdf5'%(savepath,identifier), load_parameters=2,
                load_results=0, load_derived_parameters=0, force=True)
    
    # We'll simply use auto loading so all data will be loaded when needed.
    traj.v_auto_load = True
    
    params = traj.parameters
    param1 = 'W_td_pyr'
    param2 = 'W_td_sst'
    param3 = 'W_td_vip'
    param4 = 'W_td_pv'
    
    #results = traj.res.runs.run_00000001.my_result
    # Here we load the data automatically on the fly
    SI = traj.res.summary.SI_PYR.frame
    SI_PV = traj.res.summary.SI_PV.frame
    SI_SST = traj.res.summary.SI_SST.frame
    SI_VIP = traj.res.summary.SI_VIP.frame
    
    rsc_PYRPYR = traj.res.summary.r_sc_PYRPYR.frame
    rsc_SSTSST = traj.res.summary.r_sc_SSTSST.frame
    rsc_PVPV = traj.res.summary.r_sc_PVPV.frame
    rsc_VIPVIP = traj.res.summary.r_sc_VIPVIP.frame
    rsc_PYRSST = traj.res.summary.r_sc_PYRSST.frame
    rsc_PYRPV = traj.res.summary.r_sc_PYRPV.frame
    rsc_PYRVIP = traj.res.summary.r_sc_PYRVIP.frame
    rsc_SSTPV = traj.res.summary.r_sc_SSTPV.frame
    rsc_SSTVIP = traj.res.summary.r_sc_SSTVIP.frame
    rsc_VIPPV = traj.res.summary.r_sc_VIPPV.frame
    

    
    if modeltype == 'mult':
        noTDvalue = 1.0 # 1.0 for multiplicative
        TDvalue = 2.0 # 2 for multiplicative
    elif modeltype == 'add':
        noTDvalue = 0.0 # 0.0 for additive
        TDvalue = 1.0 # 1 for additive
    
    
    rsc_dict = {
            'r_sc_PYRPYR':rsc_PYRPYR,
            'r_sc_PYRSST':rsc_PYRSST,
            'r_sc_PYRPV':rsc_PYRPV,
            'r_sc_PYRVIP':rsc_PYRVIP,
            'r_sc_SSTSST':rsc_SSTSST,
            'r_sc_SSTPV':rsc_SSTPV,
            'r_sc_SSTVIP':rsc_SSTVIP,
            'r_sc_VIPPV':rsc_VIPPV,
            'r_sc_VIPVIP':rsc_VIPVIP,
            'r_sc_PVPV':rsc_PVPV,       
            }
    SI_dict = {
            'SI_PYR':SI,
            'SI_SST':SI_SST,
            'SI_VIP':SI_VIP,
            'SI_PV':SI_PV 
    }
    
    
    
    for variable in ['rsc','SI']:
        xlabels = []
        ylabels = []
    
        if variable == 'rsc':
            data = rsc_dict
            vmin =  -.5
            vmax = .5
        elif variable == 'SI':
            data = SI_dict
            if modeltype == 'mult':
                vmin = -10
                vmax = 240
            elif modeltype == 'add':
                vmin = -10
                vmax = 100
            else:
                raise ValueError
        else: 
            raise ValueError

            
            
        mods = [[TDvalue,noTDvalue,noTDvalue,noTDvalue],[noTDvalue,TDvalue,noTDvalue,noTDvalue],[noTDvalue,noTDvalue,TDvalue,noTDvalue],[noTDvalue,noTDvalue,noTDvalue,TDvalue],[TDvalue,TDvalue,noTDvalue,noTDvalue],[TDvalue,noTDvalue,TDvalue,noTDvalue],[TDvalue,noTDvalue,noTDvalue,TDvalue],[noTDvalue,TDvalue,TDvalue,noTDvalue],[noTDvalue,TDvalue,noTDvalue,TDvalue],[noTDvalue,noTDvalue,TDvalue,TDvalue],[TDvalue,TDvalue,TDvalue,noTDvalue],[noTDvalue,TDvalue,TDvalue,TDvalue],[TDvalue,TDvalue,noTDvalue,TDvalue],[TDvalue,noTDvalue,TDvalue,TDvalue],[TDvalue,TDvalue,TDvalue,TDvalue]]
        mod_names = ['PYR','SOM','VIP','PV','PYR,SOM','PYR,VIP','PYR,PV','SOM,VIP','SOM,PV','VIP,PV','PYR,SOM,VIP','SOM,VIP,PV','PYR,SOM,PV','PYR,VIP,PV','PYR,SOM,VIP,PV']
        pairs = ['PV-PV','PYR-PV','PYR-PYR','PYR-SOM','PYR-VIP','SOM-PV','SOM-SOM','SOM-VIP','VIP-PV','VIP-VIP']
        classes = ['PV','PYR','SOM','VIP']
        change_array=np.zeros((len(mods),len(data)))
        changes = {}
        innerdict = {}
        for i, mod in enumerate(mods):
            j = 0
            xlabels.append(mod_names[i])
        
            for name, meas in sorted(data.items()):
                without_TD = meas[(meas['W_td_pyr']==noTDvalue) & (meas['W_td_sst']==noTDvalue) & (meas['W_td_vip']==noTDvalue) & (meas['W_td_pv']==noTDvalue)][name].values[0]
                with_TD = meas[(meas['W_td_pyr']==mod[0]) & (meas['W_td_sst']==mod[1]) & (meas['W_td_vip']==mod[2]) & (meas['W_td_pv']==mod[3])][name].values[0]
                change_array[i,j] = with_TD-without_TD
                changes[str(mod)] = innerdict
                innerdict[name]=with_TD-without_TD
                if i == 0:
                    ylabels.append(name)
                if i ==4:
                    if j == 1:
                        saveTD = with_TD
                        savenoTD = without_TD
                j+=1 
        
        orig_cmap = matplotlib.cm.Spectral
        mp = 1 - vmax / (vmax + abs(vmin))
        shifted_cmap = shiftedColorMap(orig_cmap, midpoint=mp, name='shifted')
        
        if plot == True:
            plt.figure(figsize=(5,4))
            plt.imshow(change_array, interpolation = 'nearest', vmin=vmin,vmax=vmax,cmap=shifted_cmap)
            
            plt.ylabel('targets of modulation')
            plt.yticks(np.arange(len(mod_names)), mod_names)
            cbar = plt.colorbar()
            
            if variable == 'rsc':
                plt.xlabel('population pair')
                plt.xticks(np.arange(len(pairs)), pairs, rotation='vertical')
                cbar.set_label('Change in noise correlation')
            
            else:
                plt.xlabel('population')
                plt.xticks(np.arange(len(classes)), classes, rotation='vertical')
                cbar.set_label('Change in selectivity')
            plt.tight_layout()
            plt.savefig('%s%s_%s.pdf'%(savepath,'%s_matrix'%variable,identifier))#, bbox_extra_artists=(lgd,), bbox_inches='tight',format='pdf', transparent=True) 
            print('saved %s%s_%s.pdf'%(savepath,'%s_matrix'%variable,identifier))
        
        changes = change_array
        scipy.io.savemat('%s%s_%s.mat'%(savepath,variable,identifier), dict(changes=changes, xlabels=xlabels, ylabels=ylabels))
        print('saved %s%s_%s.mat'%(savepath,variable,identifier))


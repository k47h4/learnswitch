#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 17:53:06 2018

plots Fig. 6B from the manuscript

@author: kwilmes
"""

from pypet import Trajectory
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps
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


    ax = var.pivot(index=other_param,columns=varied_param,values=varname).plot(colormap=cmaps.viridis)#color=cmaps.viridis(param_value/2.0))
    if var2 is not None:
        var2.pivot(index=other_param,columns=varied_param,values=var2name).plot(colormap=cmaps.magma,ax=ax)#color=cmaps.viridis(param_value/2.0))

    plt.tight_layout()
    plt.xlabel(other_param)#, fontsize=20)
    plt.ylabel(varname)#, fontsize=20)
    plt.tight_layout
    plt.savefig('%s/%s%s_with_%s_and_%s_%s.pdf'%(savepath,varname,var2name,varied_param, other_param, identifier))#, bbox_extra_artists=(lgd,), bbox_inches='tight',format='pdf', transparent=True) 

def tsplot(ax,mean,std,**kwargs):
    x = np.arange(mean.shape[0])
    cis = (mean - std, mean + std)
    ax.fill_between(x,cis[0],cis[1],alpha=0.2,**kwargs)
    ax.plot(x,mean,lw=2,**kwargs)
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

identifier = 'multiplicative'
savepath = './'

traj = Trajectory('%s'%(identifier), add_time=False)

# Let's load the trajectory from the file
traj.f_load(filename='%s/%s.hdf5'%(savepath,identifier), index = 0, load_parameters=2,
            load_results=0, load_derived_parameters=0, force=True)

traj.v_auto_load = True

params = traj.parameters
print(params)
param1 = 'W_td_pyr'
param2 = 'W_td_sst'
param3 = 'W_td_vip'
param4 = 'W_td_pv'

# Here we load the data 
rate_PYR = traj.res.summary.rate_PYR.frame
rate_PYR_std = traj.res.summary.rate_PYR_std.frame
rate_PV = traj.res.summary.rate_PV.frame
rate_PV_std = traj.res.summary.rate_PV_std.frame
rate_SST = traj.res.summary.rate_SST.frame
rate_SST_std = traj.res.summary.rate_SST_std.frame
rate_VIP = traj.res.summary.rate_VIP.frame
rate_VIP_std = traj.res.summary.rate_VIP_std.frame

peak_PYR_std = traj.res.summary.peak_PYR_std.frame
peak_PV_std = traj.res.summary.peak_PV_std.frame
peak_SST_std = traj.res.summary.peak_SST_std.frame
peak_VIP_std = traj.res.summary.peak_VIP_std.frame


noTDvalue = 1.0 
TDvalue = 1.7 

# get the means and standard deviations of the activities of the different populations
rate_PYR_woTD = rate_PYR[(rate_PYR['W_td_pyr']==noTDvalue) & (rate_PYR['W_td_sst']==noTDvalue) & (rate_PYR['W_td_vip']==noTDvalue) & (rate_PYR['W_td_pv']==noTDvalue)]['rate_PYR'].values[0]
rate_PYR_std_woTD = rate_PYR_std[(rate_PYR_std['W_td_pyr']==noTDvalue) & (rate_PYR_std['W_td_sst']==noTDvalue) & (rate_PYR_std['W_td_vip']==noTDvalue) & (rate_PYR_std['W_td_pv']==noTDvalue)]['rate_PYR_std'].values[0]
rate_PV_woTD = rate_PV[(rate_PV['W_td_pyr']==noTDvalue) & (rate_PV['W_td_sst']==noTDvalue) & (rate_PV['W_td_vip']==noTDvalue) & (rate_PV['W_td_pv']==noTDvalue)]['rate_PV'].values[0][:,0]
rate_PV_std_woTD = rate_PV_std[(rate_PV_std['W_td_pyr']==noTDvalue) & (rate_PV_std['W_td_sst']==noTDvalue) & (rate_PV_std['W_td_vip']==noTDvalue) & (rate_PV_std['W_td_pv']==noTDvalue)]['rate_PV_std'].values[0][:,0]
rate_SST_woTD = rate_SST[(rate_SST['W_td_pyr']==noTDvalue) & (rate_SST['W_td_sst']==noTDvalue) & (rate_SST['W_td_vip']==noTDvalue) & (rate_SST['W_td_pv']==noTDvalue)]['rate_SST'].values[0][:,0]
rate_SST_std_woTD = rate_SST_std[(rate_SST_std['W_td_pyr']==noTDvalue) & (rate_SST_std['W_td_sst']==noTDvalue) & (rate_SST_std['W_td_vip']==noTDvalue) & (rate_SST_std['W_td_pv']==noTDvalue)]['rate_SST_std'].values[0][:,0]
rate_VIP_woTD = rate_VIP[(rate_VIP['W_td_pyr']==noTDvalue) & (rate_VIP['W_td_sst']==noTDvalue) & (rate_VIP['W_td_vip']==noTDvalue) & (rate_VIP['W_td_pv']==noTDvalue)]['rate_VIP'].values[0][:,0]
rate_VIP_std_woTD = rate_VIP_std[(rate_VIP_std['W_td_pyr']==noTDvalue) & (rate_VIP_std['W_td_sst']==noTDvalue) & (rate_VIP_std['W_td_vip']==noTDvalue) & (rate_VIP_std['W_td_pv']==noTDvalue)]['rate_VIP_std'].values[0][:,0]

# get the peaks
peak_PYR_std_woTD = peak_PYR_std[(peak_PYR_std['W_td_pyr']==noTDvalue) & (peak_PYR_std['W_td_sst']==noTDvalue) & (peak_PYR_std['W_td_vip']==noTDvalue) & (peak_PYR_std['W_td_pv']==noTDvalue)]['peak_PYR_std'].values[0]
peak_PV_std_woTD = peak_PV_std[(peak_PV_std['W_td_pyr']==noTDvalue) & (peak_PV_std['W_td_sst']==noTDvalue) & (peak_PV_std['W_td_vip']==noTDvalue) & (peak_PV_std['W_td_pv']==noTDvalue)]['peak_PV_std'].values[0]
peak_SST_std_woTD = peak_SST_std[(peak_SST_std['W_td_pyr']==noTDvalue) & (peak_SST_std['W_td_sst']==noTDvalue) & (peak_SST_std['W_td_vip']==noTDvalue) & (peak_SST_std['W_td_pv']==noTDvalue)]['peak_SST_std'].values[0]
peak_VIP_std_woTD = peak_VIP_std[(peak_VIP_std['W_td_pyr']==noTDvalue) & (peak_VIP_std['W_td_sst']==noTDvalue) & (peak_VIP_std['W_td_vip']==noTDvalue) & (peak_VIP_std['W_td_pv']==noTDvalue)]['peak_VIP_std'].values[0]


# finally the plot of the means and s.d.s:
plt.figure(figsize=(2,2))
ax = plt.subplot(1,1,1)
tsplot(ax,rate_PYR_woTD[1100:1350],rate_PYR_std_woTD[1100:1350],label='PYR', color = 'k')
tsplot(ax,rate_PV_woTD[1100:1350],rate_PV_std_woTD[1100:1350],label='PV', color = 'orange')
tsplot(ax,rate_SST_woTD[1100:1350],rate_SST_std_woTD[1100:1350],label='SOM', color = 'b')
tsplot(ax,rate_VIP_woTD[1100:1350],rate_VIP_std_woTD[1100:1350],label='VIP', color = 'g')
lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False, framealpha=1)
plt.xlabel('time [s]')
plt.ylabel('activity')
plt.xticks(np.arange(0,250,100),np.arange(0,2.5,1))
plt.ylim(0,18)
plt.savefig('%s%s.pdf'%(savepath,identifier), bbox_extra_artists=(lgd,), bbox_inches='tight') 
print('saved to %s%s.pdf'%(savepath,identifier))

scipy.io.savemat('%s%s.mat'%(savepath,'Fig6b_responses'), {'means':{'PYR':rate_PYR_woTD[1100:1350],
                                                         'PV':rate_PV_woTD[1100:1350],
                                                         'SST':rate_SST_woTD[1100:1350],
                                                         'VIP':rate_VIP_woTD[1100:1350]},
                                                        'std':{'PYR':rate_PYR_std_woTD[1100:1350],
                                                         'PV':rate_PV_std_woTD[1100:1350],
                                                         'SST':rate_SST_std_woTD[1100:1350],
                                                         'VIP':rate_VIP_std_woTD[1100:1350]},
                                                        'time':np.arange(0,2.5,.01),
                                                        'time_unit':'seconds'
                                                         })
print('saved to %s%s.mat'%(savepath,'Fig6B_responses'))

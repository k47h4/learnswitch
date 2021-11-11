#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 17:53:06 2018

@author: kwilmes
"""

from pypet import Trajectory
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps
import seaborn as sns; sns.set(color_codes=True)
import pandas as pd
import numpy as np
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
    plt.figure(figsize=(1,1))
    ax = plt.subplot(2,1,1)
    ##Let's iterate through the columns and plot the different firing rates :
    #for param_value in SI[other_param]:
    print(var[varied_param])
    print(var[varname])
    ax = var.pivot(index=other_param,columns=varied_param,values=varname).plot(colormap=cmaps.viridis)#color=cmaps.viridis(param_value/2.0))
    if var2 is not None:
        var2.pivot(index=other_param,columns=varied_param,values=var2name).plot(colormap=cmaps.magma,ax=ax)#color=cmaps.viridis(param_value/2.0))
    #ax.set_yscale('log')
    lgd = plt.legend(title = varied_param, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False, framealpha=1)
    #plt.ylim(0,4.5)
    #lgd = plt.legend(title = "TD to SST", bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False, framealpha=1)

    plt.xlabel(other_param)
    plt.ylabel(varname)
    #plt.ylabel()

    plt.tight_layout()
    plt.savefig('%s/%s%s_with_%s_and_%s_%s.pdf'%(savepath,varname,var2name,varied_param, other_param, identifier), bbox_extra_artists=(lgd,), bbox_inches='tight',format='pdf') 


for identifier in ['pyrsst']:
    savepath = './'
    
    # This time we don't need an environment since we just going to look
    # at data in the trajectory
    traj = Trajectory('%s'%identifier, add_time=False)
    
    # Let's load the trajectory from the file
    # Only load the parameters, we will load the results on the fly as we need them
    traj.f_load(filename='%s%s.hdf5'%(savepath,identifier), load_parameters=2,
                load_results=0, load_derived_parameters=0, force = True)
    
    # We'll simply use auto loading so all data will be loaded when needed.
    traj.v_auto_load = True
    
    params = traj.parameters
    print(params)
    varied_param = 'W_td_pyr'
    other_param = 'W_td_sst'
    
    
    #results = traj.res.runs.run_00000001.my_result
    # Here we load the data automatically on the fly
    SI = traj.res.summary.SI_PYR.frame
    SI_PV = traj.res.summary.SI_PV.frame
    SI_SST = traj.res.summary.SI_SST.frame
    SI_VIP = traj.res.summary.SI_VIP.frame
    
    DIFF = traj.res.summary.DIFF_PYR.frame
    pooled_sd = traj.res.summary.pooled_sd.frame
    pooled_sd_PV = traj.res.summary.pooled_sd_PV.frame
    pooled_sd_SST = traj.res.summary.pooled_sd_SST.frame
    pooled_sd_VIP = traj.res.summary.pooled_sd_VIP.frame
    
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
    
    stdx_PYR = traj.res.summary.stdx_PYRVIP.frame
    stdx_PYRctrl = traj.res.summary.stdx_PYRPYR.frame
    stdy_VIP = traj.res.summary.stdy_PYRVIP.frame
    nc_nonnorm_PYRVIP = traj.res.summary.nc_nonnorm_PYRVIP.frame
    nc_nonnorm_PYRSST = traj.res.summary.nc_nonnorm_PYRSST.frame
    nc_nonnorm_PYRPYR = traj.res.summary.nc_nonnorm_PYRPYR.frame
    nc_nonnorm_PYRPV = traj.res.summary.nc_nonnorm_PYRPV.frame
    nc_nonnorm_VIPVIP = traj.res.summary.nc_nonnorm_VIPVIP.frame
    nc_nonnorm_SSTSST = traj.res.summary.nc_nonnorm_SSTSST.frame
    

    all_rscs = dict(PYRPYR = [rsc_PYRPYR['r_sc_PYRPYR'][0],rsc_PYRPYR['r_sc_PYRPYR'][1]],
        PYRPV = [rsc_PYRPV['r_sc_PYRPV'][0],rsc_PYRPV['r_sc_PYRPV'][1]],
        PYRSST = [rsc_PYRSST['r_sc_PYRSST'][0],rsc_PYRSST['r_sc_PYRSST'][1]],
        PYRVIP = [rsc_PYRVIP['r_sc_PYRVIP'][0],rsc_PYRVIP['r_sc_PYRVIP'][1]],        
        SSTSST = [rsc_SSTSST['r_sc_SSTSST'][0],rsc_SSTSST['r_sc_SSTSST'][1]],
        VIPVIP = [rsc_VIPVIP['r_sc_VIPVIP'][0],rsc_VIPVIP['r_sc_VIPVIP'][1]], 
        SSTPV = [rsc_SSTPV['r_sc_SSTPV'][0],rsc_SSTPV['r_sc_SSTPV'][1]],
        SSTVIP = [rsc_SSTVIP['r_sc_SSTVIP'][0],rsc_SSTVIP['r_sc_SSTVIP'][1]],
        VIPPV = [rsc_VIPPV['r_sc_VIPPV'][0],rsc_VIPPV['r_sc_VIPPV'][1]],
        PVPV = [rsc_PVPV['r_sc_PVPV'][0],rsc_PVPV['r_sc_PVPV'][1]]        
        )
    labels = ['PYR-PYR','PYR-PV','PYR-SOM','PYR-VIP','SOM-SOM','VIP-VIP','SOM-PV','SOM-VIP','VIP-PV','PV-PV']

        
    parameter = 'W_td_pyr'
    parameter2 = 'W_td_sst'
    

    conditions = ['ignore', 'attend']
    x_pos = np.arange(len(conditions))

    peaks = [SI[SI[parameter2]==1]['SI_PYR'][0],SI[SI[parameter2]==2.2]['SI_PYR'][1]] # position 0 and position 1 for the two runs without (0) and with (1) modulation
    peaksPV = [SI_PV[SI_PV[parameter2]==1]['SI_PV'][0],SI_PV[SI_PV[parameter2]==2.2]['SI_PV'][1]]
    peaksSST = [SI_SST[SI_SST[parameter2]==1]['SI_SST'][0],SI_SST[SI_SST[parameter2]==2.2]['SI_SST'][1]]
    peaksVIP = [SI_VIP[SI_VIP[parameter2]==1]['SI_VIP'][0],SI_VIP[SI_VIP[parameter2]==2.2]['SI_VIP'][1]]
    

    fig, ((ax1),(ax2),(ax3),(ax4)) = plt.subplots(1,4,figsize=(6.5,2))#, gridspec_kw = {'width_ratios':[1.2, 1]})
    ax1.bar(x_pos, peaks, width = .5, align='center', edgecolor='black', color='black', capsize=10)
    ax1.set_ylabel('SI')
    ax1.set_ylim(0,350)
    
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(conditions)
    ax1.set_title('PYR')
    #ax1.yaxis.grid(True)
    ax2.bar(x_pos, peaksPV, width = .5, align='center', edgecolor='black', color='orange', capsize=10)
    #ax2.set_ylabel('SI')
    ax2.set_ylim(0,350)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(conditions)
    ax2.set_title('PV')
    #ax2.yaxis.grid(True)
    ax3.bar(x_pos, peaksSST, width = .5, align='center', edgecolor='black', color='blue', capsize=10)
    #ax3.set_ylabel('SI')
    ax3.set_ylim(0,350)
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(conditions)
    ax3.set_title('SOM')
    #ax3.yaxis.grid(True)
    ax4.bar(x_pos, peaksVIP, width = .5, align='center', edgecolor='black', color='green', capsize=10)
    #ax3.set_ylabel('SI')
    ax4.set_ylim(0,350)
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(conditions)
    ax4.set_title('VIP')
    #ax4.yaxis.grid(True)
    plt.tight_layout()
    plt.savefig('%s/%s_%s.pdf'%(savepath,'barplot',identifier))

    all_peaks = dict(PYR=peaks, PV=peaksPV, SOM=peaksSST, VIP=peaksVIP)

    scipy.io.savemat('%sFig6D.mat'%(savepath), dict(peaks=all_peaks, xlabels=conditions))
    scipy.io.savemat('%sFig6E.mat'%(savepath), dict(rsc=all_rscs, SOM_modulation = [1.0,2.2], labels=labels))
    print('saved %sFig6D.mat and %sFig6E.mat'%(savepath,savepath))
    
    

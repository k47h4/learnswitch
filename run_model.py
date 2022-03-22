#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
L2/3 rate model consisting of pyramidal, SOM PV, and VIP populations

Created on Mon Mar  6 14:12:15 2017

@author: kwilmes
"""
import os, sys
import shutil
from tempfile import mkdtemp
import numpy as np
import pickle
#scriptpath = "/home/kwilmes/"
#sys.path.append(os.path.abspath(scriptpath))
import colormaps as cmaps
import matplotlib.cm as cmap
import seaborn as sns; sns.set(color_codes=True)
sns.set_style('darkgrid')

import pandas as pd
from OUprocess import OU_process
import random as random
import matplotlib.pyplot as plt
from sacred import Experiment
from pypet import Environment, cartesian_product
from scipy.special import expit

ex = Experiment("L23_network") 


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)
        
class TmpExpDir(object):
        """A context manager that creates and deletes temporary directories.
        """

        def __init__(self, base_dir="./"):
                self._base_dir = base_dir
                self._exp_dir = None
        
        def __enter__(self):
                # create a temporary directory into which we will store all our files
                # it will be placed into the current directory but you could change that
                self._exp_dir = mkdtemp(dir=self._base_dir)
                return self._exp_dir

        def __exit__(self, *args):
                # at the very end of the run delete the temporary directory
                # sacred will have taken care of copying all the results files over
                # to the run directoy
                if self._exp_dir is not None:
                    shutil.rmtree(self._exp_dir)



def getRandomDrift(theta=1, mu=1, sigma=1, simTime=100, steps_per_ms=100,N=1):
    OrnsteinUhlenbeck = OU_process(theta = theta, mu = mu, sigma = sigma, startPosition = mu)
    driftWalk = OrnsteinUhlenbeck.sample_path(np.linspace(0,simTime,simTime*steps_per_ms),N=N)
    return driftWalk#[0]


def get_activation(_x,r_max):
    r_0 = 1.0
    _x[_x<=0] = 0 # rectify at 0
    _x[_x>0] = (r_max-r_0)*np.tanh(_x[_x>0]/(r_max-r_0)) # tanh activation with maximum at r_max
    return _x

def calc_noise_correlations(x, y, ntrials):
    r_sc = 1/(ntrials-1) * np.sum(np.dot((x-np.mean(x))/np.std(x),(y-np.mean(y))/np.std(y)))
    r_sc_nonnorm = 1/(ntrials-1) * np.sum(np.dot((x-np.mean(x)),(y-np.mean(y))))
    return r_sc, r_sc_nonnorm, np.std(x), np.std(y)

def get_particular_stimulus_times(stimuli, no_stimuli, upto = None, startat = None):
    first = np.zeros((no_stimuli))
    endofsth = np.zeros((no_stimuli))
    startofsth = np.zeros((no_stimuli))
    last = np.zeros((no_stimuli)) 
    for i in range(no_stimuli):
        print(np.nonzero(stimuli==i)[0])
        alltimes = np.nonzero(stimuli==i)[0]
        alltimes = alltimes[alltimes>=50]
        first[i] = alltimes[0]
        if not upto == None:
            endofsth[i] = np.nonzero(stimuli[:int(upto)]==i)[0][-1]
        if startat != None:
            if startat < len(stimuli):
                startofsth[i] = np.nonzero(stimuli[int(startat):]==i)[0][0]+startat
            else:
                raise ValueError('startat is not supposed to be at the end of the simulation')
        last[i] = np.nonzero(stimuli==i)[0][-1]
    return first, endofsth, startofsth, last


def get_response(firing_rate, stimuli, no_stimuli, stim_change_time, first=None, last=None):
    if np.shape(np.shape(firing_rate))[0] == 1:
        N_neurons = 1
    else:
        N_neurons = np.shape(firing_rate)[1]
    
    #get firing rate for neurons at time of each stimulus when it is presented first (between first and first+stim_change_time) 
    #and last (between last-stim_change_time+1 and last+1):
    response = np.zeros((N_neurons,no_stimuli))
    for k in range(N_neurons):
        for i in range(no_stimuli):
            if not first is None:
                response[k,i] = np.mean(firing_rate[int(first[i]):int(first[i]+stim_change_time),k])
                print('first')
                print(int(first[i]))
                print(int(first[i]+stim_change_time))
            elif not last is None:
                response[k,i] = np.mean(firing_rate[int(last[i]-stim_change_time+1):int(last[i]+1),k])
                print('last')
                print(int(last[i]-stim_change_time+1))
                print(int(last[i]+1))
                # add +1 because: e.g. if last is 499, then stimulus goes from 450 to 499, indexing to get those rates:[450:500]  
            else:
                raise ValueError('no time given')
    return response

def get_change_in_tuning(firing_rate, stimuli, no_stimuli, stim_change_time, first, last):
    
    #get firing rate for neurons at time of each stimulus when it is presented first (between first and first+stim_change_time) 
    #and last (between last-stim_change_time+1 and last+1):
    first_response = get_response(firing_rate, stimuli, no_stimuli, stim_change_time, first=first, last=None)    
    last_response = get_response(firing_rate, stimuli, no_stimuli, stim_change_time, first=None, last=last)
    return first_response, last_response



def plot_tuning(first_resp, last_resp=None, name = ''):
    
    plt.figure()
    plt.plot([0,1,2,3],[first_resp[:,0],first_resp[:,1],first_resp[:,2],first_resp[:,3]], 'k')
    if not last_resp == None:
        plt.plot([0,1,2,3],[last_resp[:,0],last_resp[:,1],last_resp[:,2],last_resp[:,3]], 'r')
    plt.xlabel('orientation')
    plt.xlabel('firing rate')
    plt.title('%s'%name)

    plt.tight_layout()
    plt.show() 

    plt.figure()
    plt.plot([0,1,2,3,4,5,6,7],[first_resp[:,0],first_resp[:,1],first_resp[:,2],first_resp[:,3],first_resp[:,4],first_resp[:,5],first_resp[:,6],first_resp[:,7]], 'k')
    if not last_resp == None:
        plt.plot([0,1,2,3,4,5,6,7],[last_resp[:,0],last_resp[:,1],last_resp[:,2],last_resp[:,3],last_resp[:,4],last_resp[:,5],last_resp[:,6],last_resp[:,7]], 'r')
    plt.xlabel('orientation')
    plt.xlabel('firing rate')
    plt.title('%s'%name)

    plt.tight_layout()
    plt.show() 




def tsplot(ax, data,**kw):
    x = np.arange(data.shape[1])
    est = np.mean(data, axis=0)
    sd = np.std(data, axis=0)
    cis = (est - sd, est + sd)
    ax.fill_between(x,cis[0],cis[1],alpha=0.2, **kw)
    ax.plot(x,est,**kw,lw=2)
    ax.margins(x=0)

@ex.command
def run_network(params,_run):
    p = Struct(**params)
    #savepath = '/mnt/DATA/kwilmes/learningvsswitching/hdf5/'

    nsteps = int(p.simtime/p.dt)
    input_starttime = int(p.input_starttime/p.dt)
    nonpreferred_starttime = int(p.nonpreferred_starttime/p.dt)
    nonpreferred_endtime = int(p.nonpreferred_endtime/p.dt)
    modulation_starttime = int(p.modulation_starttime/p.dt)
    modulation_endtime = int(p.modulation_endtime/p.dt)
    input_endtime = int(p.input_endtime/p.dt)
    sample_resolution = int(p.sample_resolution/p.dt)

    ntrials = p.ntrials
            
    N_pyr = p.N_pyr
    N_pv = p.N_pv
    N_sst = p.N_sst
    N_vip = p.N_vip
    N_td = p.N_td
    N_ff = p.N_ff
    
    r_pyr = np.zeros((1,N_pyr))
    r_sst = np.zeros((1,N_sst))
    r_vip = np.zeros((1,N_vip))
    r_pv = np.zeros((1,N_pv))
    r_td = np.zeros((1,N_td))

    W_pyr = np.ones((N_pyr, N_pyr)) * p.W_pyr 
    W_sst_pyr = np.ones((N_sst, N_pyr)) * p.W_sst_pyr    
    W_pv_pyr = np.ones((N_pv, N_pyr)) * p.W_pv_pyr 
    W_pyr_sst = np.ones((N_pyr, N_sst)) * p.W_pyr_sst    
    
    W_vip_sst = np.ones((N_vip, N_sst)) * p.W_vip_sst 
    W_vip_pv = np.ones((N_vip, N_pv)) * p.W_vip_pv 
    W_vip_pyr = np.ones((N_vip, N_pyr)) * p.W_vip_pyr 
    
    W_pyr_vip = np.ones((N_pyr, N_vip)) * p.W_pyr_vip
    W_sst_vip = np.ones((N_sst, N_vip)) * p.W_sst_vip 

    W_pyr_pv = np.ones((N_pyr, N_pv)) * p.W_pyr_pv    
    W_sst_pv = np.ones((N_sst, N_pv)) * p.W_sst_pv  
    W_pv_pv = np.ones((N_pv, N_pv)) * p.W_pv_pv 
    W_vip_vip = np.ones((N_vip, N_vip)) * p.W_vip_vip  
    W_pv_vip = np.ones((N_pv, N_vip)) * p.W_pv_vip  
    
    W_td_pyr = p.W_td_pyr
    W_td_vip = p.W_td_vip
    W_td_sst = p.W_td_sst
    W_td_pv = p.W_td_pv


    #print('td modulation of pyr, sst, vip')
    #print(W_td_pyr)
    #print(W_td_sst)
    #print(W_td_vip)
    W_ff_pyr = p.W_ff_pyr
    W_ff_sst = p.W_ff_sst
    W_ff_pv = p.W_ff_pv
    ff = p.ff
    td = p.td
    xi_ff = p.xi_ff
    xi_td = p.xi_ff
    xi_td_pyr = 0.0
    xi_td_pv = 0.0
    xi_td_sst = 0.0
    xi_td_vip = 0.0
    
    max_rate = p.max_rate
    threshold = p.threshold
    dispersion = p.dispersion
    

    rates_pyr = np.zeros((int(nsteps/sample_resolution),N_pyr,ntrials))
    rates_sst = np.zeros((int(nsteps/sample_resolution),N_sst,ntrials))
    rates_vip = np.zeros((int(nsteps/sample_resolution),N_vip,ntrials))
    rates_pv = np.zeros((int(nsteps/sample_resolution),N_pv,ntrials))
    rates_ff = np.zeros((int(nsteps/sample_resolution),ntrials))
    rates_td = np.zeros((int(nsteps/sample_resolution),ntrials))

    
    for trial in range(ntrials):
                
        np.random.seed(trial)

        for t in range(nsteps):
    
            if t > nonpreferred_starttime and t < nonpreferred_endtime:
                ff_noise = np.random.randn(1)
                r_ff = 0.2
            elif t > input_starttime and t < input_endtime: # during feedforward input
                ff_noise = np.random.randn(1)
                r_ff = 1.0                    
            else:
                r_ff = 0
                ff_noise = 0
    
            if ((t > modulation_starttime) and (t < modulation_endtime)): # during modulation time
                r_td = 1.0
                td_noise= np.random.randn(1)
                if p.modulation == 'additive':                    
                    sigma_pyr = p.sigma 
                    sigma_pv = p.sigma 
                    sigma_sst = p.sigma 
                    sigma_vip = p.sigma 
                    xi_td_pyr = W_td_pyr * xi_td
                    xi_td_pv = W_td_pv * xi_td
                    xi_td_sst = W_td_sst * xi_td
                    xi_td_vip = W_td_vip * xi_td                    
                elif p.modulation == 'multiplicative':
                    sigma_pyr = p.sigma #
                    sigma_pv = p.sigma #
                    sigma_sst = p.sigma #
                    sigma_vip = p.sigma #
                    xi_td_pyr = (W_td_pyr-1) * xi_td
                    xi_td_pv = (W_td_pv-1) * xi_td
                    xi_td_sst = (W_td_sst-1) * xi_td
                    xi_td_vip = (W_td_vip-1) * xi_td                                                     
                else:
                    sigma_pyr = p.sigma #
                    sigma_pv = p.sigma #
                    sigma_sst = p.sigma #
                    sigma_vip = p.sigma #
                    xi_td_pyr = 0
                    xi_td_pv = 0
                    xi_td_sst = 0
                    xi_td_vip = 0    
            else:
                r_td = 0.0
                td_noise = 0
                sigma_pyr = p.sigma
                sigma_pv = p.sigma
                sigma_sst = p.sigma
                sigma_vip = p.sigma
                
                
            k = 1
            n = 1
            if p.modulation == 'additive':            
                r_pyr += p.dt*(-1*r_pyr + get_activation(np.dot(r_pyr,W_pyr) - np.dot(r_pv,W_pv_pyr) - np.dot(r_sst,W_sst_pyr) - np.dot(r_vip,W_vip_pyr) + np.dot(r_ff,W_ff_pyr) + np.dot(r_td,W_td_pyr) + p.baseline_pyr + (np.sqrt(p.dt)/p.dt) * (sigma_pyr * ((1-xi_ff-xi_td_pyr) * np.random.randn(2) + xi_ff * ff_noise + (xi_td_pyr) * td_noise)), p.r_max))/p.tau_E # add noise: np.random.randn(1) * p.sigma    
                r_pv += p.dt*(-1*r_pv + get_activation(np.dot(r_pyr,W_pyr_pv) - np.dot(r_pv,W_pv_pv) - np.dot(r_sst,W_sst_pv) - np.dot(r_vip,W_vip_pv) + np.dot(r_ff,W_ff_pv) + np.dot(r_td,W_td_pv) + p.baseline_pv + (np.sqrt(p.dt)/p.dt) * (sigma_pv * ((1-xi_ff-xi_td_pv) * np.random.randn(2) + xi_ff * ff_noise + xi_td_pv * td_noise)), p.r_max_i))/p.tau_I # add noise: np.random.randn(1) * p.sigma 
                r_sst += p.dt*(-1*r_sst + get_activation(np.dot(r_pyr,W_pyr_sst) - np.dot(r_vip,W_vip_sst) + np.dot(r_td,W_td_sst) + np.dot(r_ff,W_ff_sst) + p.baseline_sst + (np.sqrt(p.dt)/p.dt) * (sigma_sst * ((1-xi_td_sst) * np.random.randn(2) + xi_td_sst * td_noise)), p.r_max_i))/p.tau_I # add noise: np.random.randn(1) * p.sigma 
                r_vip += p.dt*(-1*r_vip + get_activation(np.dot(r_pyr,W_pyr_vip) - np.dot(r_sst,W_sst_vip) - np.dot(r_pv,W_pv_vip) - np.dot(r_pv,W_pv_vip) + np.dot(r_td,W_td_vip) + p.baseline_vip + (np.sqrt(p.dt)/p.dt) * (sigma_vip * ((1-xi_td_vip) * np.random.randn(2) + (xi_td_vip) * td_noise)), p.r_max_i))/p.tau_I # add noise: np.random.randn(1) * p.sigma 
            elif p.modulation == 'multiplicative':
                r_pyr += p.dt*(-1*r_pyr + get_activation(np.dot(r_td,W_td_pyr) * (np.dot(r_ff,W_ff_pyr) + p.baseline_pyr) + np.dot(r_pyr,W_pyr) - np.dot(r_pv,W_pv_pyr) - np.dot(r_sst,W_sst_pyr) - np.dot(r_vip,W_vip_pyr) + (np.sqrt(p.dt)/p.dt) * (sigma_pyr * ((1-xi_ff-xi_td_pyr) * np.random.randn(2) + xi_ff * ff_noise + (xi_td_pyr) * td_noise)), p.r_max))/p.tau_E # add noise: np.random.randn(1) * p.sigma    
                r_pv += p.dt*(-1*r_pv + get_activation(np.dot(r_td,W_td_pyr) * (np.dot(r_ff,W_ff_pv) + p.baseline_pv) + np.dot(r_pyr,W_pyr_pv) - np.dot(r_pv,W_pv_pv) - np.dot(r_sst,W_sst_pv) - np.dot(r_vip,W_vip_pv) + (np.sqrt(p.dt)/p.dt) * (sigma_pv * ((1-xi_ff-xi_td_pv) * np.random.randn(2) + xi_ff * ff_noise + xi_td_pv * td_noise)), p.r_max_i))/p.tau_pv # add noise: np.random.randn(1) * p.sigma 
                r_sst += p.dt*(-1*r_sst + get_activation(np.dot(r_td,W_td_sst) * p.baseline_sst + np.dot(r_pyr,W_pyr_sst) - np.dot(r_vip,W_vip_sst) + (np.sqrt(p.dt)/p.dt) * (sigma_sst * ((1-xi_td_sst) * np.random.randn(2) + xi_td_sst * td_noise)), p.r_max_i))/p.tau_sst # add noise: np.random.randn(1) * p.sigma 
                r_vip += p.dt*(-1*r_vip + get_activation(np.dot(r_td,W_td_vip) * p.baseline_vip + np.dot(r_pyr,W_pyr_vip) - np.dot(r_sst,W_sst_vip) - np.dot(r_pv,W_pv_vip) - np.dot(r_vip,W_vip_vip) + (np.sqrt(p.dt)/p.dt) * (sigma_vip * ((1-xi_td_vip) * np.random.randn(2) + (xi_td_vip) * td_noise)), p.r_max_i))/p.tau_vip # add noise: np.random.randn(1) * p.sigma 
            else:
                raise ValueError
                

            r_pyr[r_pyr<0] = 0
            r_sst[r_sst<0] = 0
            r_pv[r_pv<0] = 0
            r_vip[r_vip<0] = 0

    
            # record data:
            if t % sample_resolution == 0:
                rates_ff[int(t/sample_resolution),trial] = r_ff
                rates_td[int(t/sample_resolution),trial] = r_td
                rates_pyr[int(t/sample_resolution),:,trial] = r_pyr
                rates_sst[int(t/sample_resolution),:,trial] = r_sst
                rates_vip[int(t/sample_resolution),:,trial] = r_vip
                rates_pv[int(t/sample_resolution),:,trial] = r_pv

    rates = {
            'PYR':rates_pyr,
            'PV':rates_pv,
            'SST':rates_sst,
            'VIP':rates_vip
            }

    # start and end of preferred stimulus
    start = int(input_starttime*p.dt+0)
    end = int(input_starttime*p.dt+100)
    
    print('start')
    print(start)
    print(end)
    print(np.mean(rates['PYR'][:,0,:],1))
    # start and end of non-preferred stimulus
    np_start = int(nonpreferred_starttime*p.dt)
    np_end = int(nonpreferred_starttime*p.dt+100)

    # rates and selectivity
    responses = []
    rates_pref = {}
    for pop_name, pop_rate in rates.items():
        for i in range(2):

            rate_pref = np.mean(pop_rate[start:end,i,:],0)    
            rate_nonpref = np.mean(pop_rate[np_start:np_end,i,:],0)    
            pooled_sd = np.sqrt(((ntrials-1)*np.std(rate_pref)**2 + (ntrials-1)*np.std(rate_nonpref)**2)/(2*(ntrials-1))) 
            SI = (np.mean(rate_pref)-np.mean(rate_nonpref))/pooled_sd
            DIFF = (np.mean(rate_pref)-np.mean(rate_nonpref))
            responses.append({'population':pop_name+str(i),'rate_pref':np.mean(rate_pref),'rate_nonpref':np.mean(rate_nonpref), 'SI':SI, 'DIFF':DIFF, 'pooled_sd':pooled_sd})
            rates_pref[pop_name+str(i)] = rate_pref
    responses = pd.DataFrame(responses)
    nc = {}
    nc_nonnorm = {}
    stdx = {}
    stdy = {}

    for pop_name1, rate_pref1 in rates_pref.items():        
        for pop_name2, rate_pref2 in rates_pref.items():        
            nc_tmp, nc_nonnormtmp, stdx_tmp, stdy_tmp = calc_noise_correlations(rate_pref1, rate_pref2, ntrials)
            #print(nc_tmp)
            nc[pop_name1+pop_name2]=nc_tmp
            nc_nonnorm[pop_name1+pop_name2]=nc_nonnormtmp
            stdx[pop_name1+pop_name2]= stdx_tmp
            stdy[pop_name1+pop_name2]= stdy_tmp
    
    results = {
        'rate_PYR':np.mean(rates['PYR'][:,0,:],1),
        'rate_SST':np.mean(rates['SST'],2),
        'rate_PV':np.mean(rates['PV'],2),
        'rate_VIP':np.mean(rates['VIP'],2),
        'rate_PYR_std':np.std(rates['PYR'][:,0,:],1),
        'rate_SST_std':np.std(rates['SST'],2),
        'rate_PV_std':np.std(rates['PV'],2),
        'rate_VIP_std':np.std(rates['VIP'],2),
        'peak_PYR_std':np.std(np.max(rates['PYR'][:,0,:],0)),
        'peak_SST_std':np.std(np.max(rates['SST'],0)),
        'peak_PV_std':np.std(np.max(rates['PV'],0)),
        'peak_VIP_std':np.std(np.max(rates['VIP'],0)),
        'SI_PYR':responses[responses.population=='PYR0']['SI'].values[0],
        'DIFF_PYR':responses[responses.population=='PYR0']['DIFF'].values[0],
        'pooled_sd':responses[responses.population=='PYR0']['pooled_sd'].values[0],   
        'pooled_sd_PV':responses[responses.population=='PV0']['pooled_sd'].values[0],        
        'pooled_sd_SST':responses[responses.population=='SST0']['pooled_sd'].values[0],        
        'pooled_sd_VIP':responses[responses.population=='VIP0']['pooled_sd'].values[0],        
        'SI_SST':responses[responses.population=='SST0']['SI'].values[0],
        'SI_PV':responses[responses.population=='PV0']['SI'].values[0],
        'SI_VIP':responses[responses.population=='VIP0']['SI'].values[0],
        'pref_PYR':responses[responses.population=='PYR0']['rate_pref'].values[0],
        'pref_PV':responses[responses.population=='PV0']['rate_pref'].values[0],
        'pref_SST':responses[responses.population=='SST0']['rate_pref'].values[0],
        'pref_VIP':responses[responses.population=='VIP0']['rate_pref'].values[0],
        'non_pref_PYR':responses[responses.population=='PYR0']['rate_nonpref'].values[0],
        'non_pref_PV':responses[responses.population=='PV0']['rate_nonpref'].values[0],
        'non_pref_SST':responses[responses.population=='SST0']['rate_nonpref'].values[0],
        'non_pref_VIP':responses[responses.population=='VIP0']['rate_nonpref'].values[0],
        'r_sc_PYRPYR':nc['PYR0PYR1'],
        'r_sc_SSTSST':nc['SST0SST1'],
        'r_sc_VIPVIP':nc['VIP0VIP1'],
        'r_sc_PVPV':nc['PV0PV1'],
        'r_sc_PYRPV':nc['PYR0PV0'],
        'r_sc_PYRSST':nc['PYR0SST0'],
        'r_sc_PYRVIP':nc['PYR0VIP0'],
        'r_sc_SSTPV':nc['SST0PV0'],
        'r_sc_SSTVIP':nc['SST0VIP0'],
        'r_sc_VIPPV':nc['VIP0PV0'],
        'nc_nonnorm_PYRPYR':nc_nonnorm['PYR0PYR0'],
        'nc_nonnorm_PYRPV':nc_nonnorm['PYR0PV0'],
        'nc_nonnorm_PYRSST':nc_nonnorm['PYR0SST0'],
        'nc_nonnorm_SSTSST':nc_nonnorm['SST0SST0'],
        'nc_nonnorm_PYRVIP':nc_nonnorm['PYR0VIP0'],
        'nc_nonnorm_VIPVIP':nc_nonnorm['VIP0VIP0'],
        'stdx_PYRVIP':stdx['PYR0VIP0'],
        'stdx_PYRPYR':stdx['PYR0PYR0'],
        'stdy_PYRVIP':stdy['PYR0VIP0']
    }

 
    return results

def my_pypet_wrapper(traj, varied_params):
    
    N_pyr = 2
        
    params = {
        # simulation parameters
        'allplots':True,
        'cluster': False,
        'seed' : 5008,
        'ntrials' : 100,
        'dt' : 0.1, # time step in milliseconds
        'simtime' :1500,
        'nonpreferred_starttime':500,
        'nonpreferred_endtime':800,
        'input_starttime' : 1200,
        'input_endtime' : 1500,
        'modulation_starttime' :0,
        'modulation_endtime' :1500,
        #change this to switch between additive and multiplicative:
        'modulation': 'additive',#'multiplicative',#'additive',# options: 'additive' or 'multiplicative'
        'dispersion':1.0,
        
    
        # number of neurons
        'N_pyr' : N_pyr, # Number of excitatory L23 PYR cells
        'N_vip' : 2,     # Number of inhibitory VIP cells
        'N_sst' : 2,     # Number of inhibitory SST cells
        'N_pv' : 2,      # Number of inhibitory PV cells 
        'N_td' : 1,      # Top-down modulation
        'N_ff' : 1,      # Feedforward input
        
        #time constants
        'tau_I' : 40,
        'tau_sst' : 40,
        'tau_pv' : 40,
        'tau_vip' : 40,
        'tau_E' : 80,
        
        #maximum firing rates
        'r_max' : 20.0,#
        'r_max_i': 20.0,#
        'max_rate':10.0,# 
        'threshold':5,#
        
        # input
        'ff': 1.0,
        'td': 1.0,

        # noise
        'sigma' : .5 * np.sqrt(2),
        'xi_td': (1/3.0),
        'xi_ff': (1/3.0),        
        
        # connectivity
        'condition' : 0,
        'W_pyr' : .017, # between pyramidal cells
        'W_pyr_pv' : .8535, 
        'W_pyr_sst' : 1.285, 
        'W_pyr_vip' : 2.104,
        'W_pv_pyr' : .956,
        'W_sst_pyr' : .512,
        'W_sst_vip' : .9,
        'W_sst_pv' : .307,
        'W_pv_pv' : .99,        
        'W_pv_vip' : .184,        
        'W_vip_pyr' : .045,
        'W_vip_sst' : .14,
        'W_vip_pv' : .09,
        'W_vip_vip' : .0,

        
        'W_td_pyr' : 1.0,        
        'W_td_vip' : 1.0,        
        'W_td_sst' : 1.0,        
        'W_td_pv' : 1.0,        

        'W_ff_pyr' : np.array([[17.8,17.8]]),
        'W_ff_sst' : np.array([[0.0]]),
        'W_ff_pv' : np.array([[10.0]]),
        'baseline_pyr' : 6.0,
        'baseline_pv' : 4.0,
        'baseline_sst' : 1.2,
        'baseline_vip' : 4.6,
        
        # recording
        'sample_resolution' : 1,
        'plot':True, 
        }
    
    
    for varied_param in varied_params:
        params[varied_param] = traj[varied_param]
    
    results = run_network(params)
    for key, value in results.items():
        traj.f_add_result('%s.$'%key, value, comment='%s `run_network`'%key)
    
    return results
    
    
def postproc(traj, result_list):
    """Postprocessing, sorts results into a data frame.

    :param traj:

        Container for results and parameters

    :param result_list:

        List of tuples, where first entry is the run index and second is the actual
        result of the corresponding run.

    :return:
    """


    param1 = 'W_td_pyr'
    param2 = 'W_td_sst'
    param3 = 'W_td_vip'
    param4 = 'W_td_pv'
    param5 = 'relative_modulation'

    
    param1_range = traj.par.f_get(param1).f_get_range()
    param2_range = traj.par.f_get(param2).f_get_range()
    param3_range = traj.par.f_get(param3).f_get_range()
    param4_range = traj.par.f_get(param4).f_get_range()

    param1_index = sorted(set(param1_range))
    param2_index = sorted(set(param2_range))
    param3_index = sorted(set(param3_range))
    param4_index = sorted(set(param4_range))
    
    
    for key in result_list[0][1]:
        # create a result_frame for each measurable
        results = []
        # go trough all simulations and store the resuls in the data frame
        for result_tuple in result_list:
            run_idx = result_tuple[0]
            var = result_tuple[1][key]
            param1_val = param1_range[run_idx]
            param2_val = param2_range[run_idx]
            param3_val = param3_range[run_idx]
            param4_val = param4_range[run_idx]

            results.append({param1:param1_val, param2:param2_val, param3:param3_val, param4:param4_val, key:var}) # Put the data into the
            # data frame
        results_frame = pd.DataFrame(results)

        # Finally we going to store our results_frame into the trajectory
        #print(results_frame)

        # Finally we going to store our results_frame into the trajectory
        #print(results_frame)
        traj.f_add_result('summary.%s'%key, frame = results_frame,
                          comment='Contains a pandas data frame with all %s'%key)


@ex.automain
def main():
    #change this to switch between additive and multiplicative#
    identifier = 'additive'#'multiplicative'# options:'additive' or 'multiplicative' or 'pyrsst'
    
    savepath = './'
    if not os.path.exists(savepath):
        os.mkdir(savepath)

    # Create the environment
    env = Environment(trajectory='%si'%identifier,
                  comment='Experiment to measure selectivity '
                        'and noise correlations. '
                        ' of different cell populations'
                        'Exploring different top-down inputs, '
                        'as well as connectivity',
                  add_time=False, # We don't want to add the current time to the name,
                  log_config='DEFAULT',
                  multiproc=True,
                  ncores=1, # run code in parallel on ncores cores
                  filename=savepath, # We only pass a folder here, so the name is chosen
                  # automatically to be the same as the Trajectory)
                  )
    traj = env.traj

    # Now add the parameters and some exploration
    param1 = 'W_td_pyr'
    param2 = 'W_td_sst'
    param3 = 'W_td_vip'
    param4 = 'W_td_pv'

    varied_params = [param1,param2,param3,param4]
    
    if identifier == 'multiplicative':
 
        traj.f_add_parameter(param1, 1.0)
        traj.f_add_parameter(param2, 1.0)
        traj.f_add_parameter(param3, 1.0)
        traj.f_add_parameter(param4, 1.0)

        explore_dict = {param1: np.arange(1.0,2.1,1.0).tolist(),
                    param2: np.arange(1.0,2.1,1.0).tolist(),
                    param3: np.arange(1.0,2.1,1.0).tolist(),
                    param4: np.arange(1.0,2.1,1.0).tolist()}

        explore_dict = cartesian_product(explore_dict, (param1, param2, param3, param4))
        # The second argument, the tuple, specifies the order of the cartesian product,
        # The variable on the right most side changes fastest and defines the
        # 'inner for-loop' of the cartesian product


    elif identifier == 'additive':

        traj.f_add_parameter(param1, 0.0)
        traj.f_add_parameter(param2, 0.0)
        traj.f_add_parameter(param3, 0.0)
        traj.f_add_parameter(param4, 0.0)


        explore_dict = {param1: np.arange(0.0,1.1,1.0).tolist(),
                    param2: np.arange(0.0,1.1,1.0).tolist(),
                    param3: np.arange(0.0,1.1,1.0).tolist(),
                    param4: np.arange(0.0,1.1,1.0).tolist()}

        explore_dict = cartesian_product(explore_dict, (param1, param2, param3, param4))


    elif identifier == 'pyrsst':


        modulation_range = np.arange(0.0,1.3,1.2)
        param2_range = 1.0 + modulation_range    
        param1_range = 1.0 + 0.7 * modulation_range # PYR modulation is 0.7 * SOM modulation

        traj.f_add_parameter(param1, 1.0)
        traj.f_add_parameter(param2, 1.0) 
        traj.f_add_parameter(param3, 1.0)
        traj.f_add_parameter(param4, 1.0)

        explore_dict = {param2: param2_range.tolist(),
                        param1: param1_range.tolist(),
                        param3: np.array([1.0,1.0]).tolist(),
                        param4: np.array([1.0,1.0]).tolist()}

    else:
        raise ValueError 

    traj.f_explore(explore_dict)
    
    
    #traj.f_explore( {param1: np.arange(0.0, 2.2, 1.0).tolist() } )

    
    # Ad the postprocessing function
    env.add_postprocessing(postproc)

    # Run your wrapping function instead of your simulator
    env.run(my_pypet_wrapper,varied_params)

    #run_network()
    
    # Finally disable logging and close all log-files
    env.disable_logging()
    
    

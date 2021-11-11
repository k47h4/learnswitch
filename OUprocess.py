#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PyProcess
@author: Cameron Davidson-Pilon
"""
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.special import gammainc
import warnings
import matplotlib.pyplot as plt
import pdb

np.seterr( divide="raise")

class Diffusion_process(object):
    #
    # Class that can be overwritten in the subclasses:
    # _sample_position(t, N=1)
    # _mean(t)
    # _var(t)
    # _sample_path(times, N=1)
    #
    # Class that should be present in subclasses:
    #
    # _transition_pdf(x,t,y)
    #
    #
    
    def __init__(self, startTime = 0, startPosition = 0, endTime = None, endPosition = None):
    
        self.conditional=False
        self.startTime = startTime
        self.startPosition = startPosition
        if ( endTime != None ) and ( endPosition != None ):
            assert endTime > startTime, "invalid parameter: endTime > startTime"
            self.endTime = endTime
            self.endPosition = endPosition
            self.conditional = True
        elif ( endTime != None ) != ( endPosition != None ):
            raise Exception( "invalid parameter:", "Must include both endTime AND endPosition or neither" )

    
    def transition_pdf(self,t,y):
        self._check_time(t)
        "this method calls self._transition_pdf(x,t,y) in the subclass"
        try:
            if not self.conditional:
                return self._transition_pdf(self.startPosition, t-self.startTime, y)
            else:
                return self._transition_pdf(self.startPosition, t-self.startTime, y)*self._transition_pdf(y, self.endTime-t, self.endPosition)\
                        /self._transition_pdf(self.startPosition,self.endTime - self.startTime, self.endPosition)
        except AttributeError:
            raise AttributeError("Attn: transition density for process is not defined.")

    def expected_value(self,t, f= lambda x:x, N=1e6):
        """
This function calculates the expected value of E[ X_t | F ] where F includes start conditions and possibly end conditions.
"""
        warnings.warn( "Performing Monte Carlo with %d simulations."%N)
        self._check_time(t)
        if not self.conditional:
            return f( self.sample_position(t, N) ).mean()
        else:
            #This uses a change of measure technique.
            """
sum=0
self.conditional=False
for i in range(N):
X = self.generate_position_at(t)
sum+=self._transition_pdf(X,self.endTime-t,self.endPosition)*f(X)
self.conditional=True
"""
            x = self.sample_position(t, N)
            mean = (self._transition_pdf(x,self.endTime-t,self.endPosition)*f(x)).mean()
            return mean/self._transition_pdf(self.startPosition, self.endTime-self.startTime, self.endPosition)
        
    def sample_position(self,t, N=1):
        """
if _get_position_at() is not overwritten in a subclass, this function will use euler scheme
"""

        self._check_time(t)
        return self._sample_position(t, N)
        
    
    def mean(self,t):
        self._check_time(t)
        return self._mean(t)
    
    
    def var(self,t):
        self._check_time(t)
        return self._var(t)
    
    def sample_path(self,times, N = 1):
        self._check_time(times[0])
        return self._sample_path(times, N)
    
    
    
    
    
    def _sample_path(self,times, N=1 ):
        return self.Euler_scheme(times, N)
    
    def _var(self,t, N = 1e6):
        """
var = SampleVarStat()
for i in range(10000):
var.push(self.generate_position_at(t))
return var.get_variance()
"""
        return self.sample_position(t, n).var()
    
    def _mean(self,t):
        return self.expected_value( t)
    
    def _sample_position(self,t):
        return self.Euler_scheme(t)
        
    def _transition_pdf(self,x,t,y):
        warning.warn( "Attn: transition pdf not defined" )
    
    def _check_time(self,t):
        if t<self.startTime:
            warnings.warn( "Attn: inputed time not valid (check if beginning time is less than startTime)." )
    
    def Euler_scheme(self, times,delta=0.001):
        """
times is an array of floats.
The process needs the methods drift() and diffusion() defined.
"""
        warnings.warn("Attn: starting an Euler scheme to approxmiate process.")
        
        Nor = stats.norm()
        finalTime = times[-1]
        steps = int(finalTime/delta)
        t = self.startTime
        x=self.startPosition
        path=[]
        j=0
        time = times[j]
        for i in xrange(steps):
            if t+delta>time>t:
                delta = time-t
                x += drift(x,t)*delta + np.sqrt(delta)*diffusion(x,t)*Nor.rvs()
                path.append((x,time))
                delta=0.001
                j+=1
                time = times[j]
            else:
                x += drift(x,t)*delta + np.sqrt(delta)*diffusion(x,t)*Nor.rvs()
                t += delta
            
        return path

    def Milstein_Scheme(self, times, delta = 0.01 ):
        if not all( map( lambda x: hasattr(self, x), ( 'drift', 'diffusion', 'diffusion_prime' ) ) ):
            raise AttributeError("The process does not have 'drift', 'diffusion', or 'diffusion_prime' methods")
        
        pass
    
    
        
        
    def plot(self, times ,N=1, **kwargs):
        assert N >= 1, "N must be greater than 0."
        try:
            self._check_time(times[-1] )
            plt.plot(times, self.sample_path(times, N).T, **kwargs )
        except:
            self._check_time(times)
            times = np.linspace(self.startTime, times, 100)
            path = self.sample_path(times, N).T
            plt.plot(times, path, **kwargs )
        plt.xlabel( "time, $t$")
        plt.ylabel( "position of process")
        plt.show()
        return
        
class Wiener_process(Diffusion_process):
    """
This implements the famous Wiener process. I choose not to call it Brownian motion, as
brownian motion is a special case of this with 0 drift and variance equal to t.
dW_t = mu*dt + sigma*dB_t
W_t ~ N(mu*t, sigma**2t)
parameters:
mu: the constant drift term, float
sigma: the constant volatility term, float > 0
"""
    
    def __init__(self, mu, sigma, startTime = 0, startPosition = 0, endPosition = None, endTime = None):
        super(Wiener_process,self).__init__(startTime, startPosition, endTime, endPosition)
        self.mu = mu
        self.sigma = sigma
        self.Nor = stats.norm()
    
    def _transition_pdf(self,x,t,y):
        return np.exp(-(y-x-self.mu*(t-self.startTime))**2/(2*self.sigma**2*(t-self.startTime)))\
            /np.sqrt(2*pi*self.sigma*(t-self.startTime))
    
    def _mean(self,t):
        if self.conditional:
            delta1 = t - self.startTime
            delta2 = self.endTime - self.startTime
            return self.startPosition + self.mu*delta1 + (self.endPosition-self.startPosition-self.mu*delta2)*delta1/delta2
        else:
            return self.startPosition+self.mu*(t-self.startTime)

    def _var(self,t):
        if self.conditional:
            delta1 = self.sigma**2*(t-self.startTime)*(self.endTime-t)
            delta2 = self.endTime-self.startTime
            return delta1/delta2
        else:
            return self.sigma**2*(t-self.startTime)
        
    def _sample_position(self,t, n=1):
        """
This incorporates both conditional and unconditional
"""
        return self.mean(t) + np.sqrt(self.var(t))*self.Nor.rvs(n)

    def _sample_path(self,times, N = 1):

        path=np.zeros( (N,len(times)) )
        path[ :, 0] = self.startPosition
        times = np.insert( times, 0,self.startTime)
        deltas = np.diff( times )

        if not self.conditional:
            path += np.random.randn( N, len(times)-1 )*self.sigma*np.sqrt(deltas) + self.mu*deltas
            return path.cumsum(axis=1)
        
                
        else:
            """
Alternatively, this can be accomplished by sampling directly from a multivariate normal given a linear
projection. Ie
N | N dot 1 = endPosition ~ Nor( 0, Sigma ), where Sigma is a diagonal matrix with elements proportional to
the delta. This only problem with this is Sigma is too large for very large len(times).
"""
            T = self.endTime - self.startTime
            x = self.startTime
            for i, delta in enumerate( deltas ):
                x = x*(1-delta/T)+self.endPosition*delta/T + self.sigma*np.sqrt(delta/T*(T-delta))*self.Nor.rvs(N)
                T = T - delta
                path[:,i] = x
            if abs(T -0)<1e-10:
                path[:,-1] = self.endPosition
        return path
    
    def generate_max(self,t):
        pass
    
    def generate_min(self,t):
        pass
    
    def drift(t,x):
        return self.mu
    def diffusion(t,x):
        return self.sigma
    def diffusion_prime(t,x):
        return 0

        
        
        
#have a Vasicek model that has an alternative parameterizatiom but essentially just maps to OU_process

class OU_process(Diffusion_process):
    
    """
The Orstein-Uhlenbeck process is a mean reverting model that is often used in finance to model interest rates.
It is also known as a Vasicek model. It is defined by the SDE:
dOU_t = theta*(mu-OU_t)*dt + sigma*dB_t$
The model flucuates around its long term mean mu. mu is also a good starting Position of the process.
There exists a solution to the SDE, see wikipedia.
parameters:
theta: float > 0
mu: float
sigma: float > 0
"""
    def __init__(self, theta, mu, sigma, startTime = 0, startPosition = 0, endPosition = None, endTime = None ):
        assert sigma > 0 and theta > 0, "theta > 0 and sigma > 0."
        super(OU_process, self).__init__(startTime, startPosition, endTime, endPosition)
        self.theta = theta
        self.mu = mu
        self.sigma = sigma
        self.Normal = stats.norm()


    def _mean(self,t):
        if self.conditional:
            return super(OU_process,self)._mean(t) #TODO
        else:
            return self.startPosition*np.exp(-self.theta*(t-self.startTime))+self.mu*(1-np.exp(-self.theta*(t-self.startTime)))
                                                                
    def _var(self,t):
        if self.conditional:
            return super(OU_process,self)._get_variance_at(t)
        else:
            return self.sigma**2*(1-np.exp(-2*self.theta*t))/(2*self.theta)
    
    def _transition_pdf(self,x,t,y):
            mu = x*np.exp(-self.theta*t)+self.mu*(1-np.exp(-self.theta*t))
            sigmaSq = self.sigma**2*(1-np.exp(-self.theta*2*t))/(2*self.theta)
            return np.exp(-(y-mu)**2/(2*sigmaSq))/np.sqrt(2*pi*sigmaSq)
    
    
    def _sample_position(self,t):
        if not self.conditional:
            return self.get_mean_at(t)+np.sqrt(self.get_variance_at(t))*self.Normal.rvs()
        else:
            #this needs to be completed
            return super(OU_process,self)._generate_position_at(t)
    
    def sample_path(self,times, N= 1, return_normals = False):
        "the parameter Normals = 0 is used for the Integrated OU Process"
        if not self.conditional:
            path=np.zeros( (N,len(times)) )
            times = np.insert( times, 0,self.startTime)
            path[:,0] = self.startPosition
            
            deltas = np.diff(times)
            normals = np.random.randn( N, len(times)-1 )
            x = self.startPosition*np.ones( N)
            print(self.sigma)
            print(self.theta)
            print(deltas)
            sigma = np.sqrt(self.sigma**2*(1-np.exp(-2*self.theta*deltas))/(2*self.theta))
            for i, delta in enumerate(deltas):
                mu = self.mu + np.exp(-self.theta*delta)*(x-self.mu)
                path[:, i] = mu + sigma[i]*normals[:,i]
                x = path[:,i]
                """
It would be really cool if there was a numpy func like np.cumf( func, array )
that applies func(next_x, prev_x) to each element. For example, lambda x,y: y + x
is the cumsum, and lambda x,y: x*y is the cumprod function.
"""
            if return_normals:
                return (path, normals )
            else:
                return path
        
        else:
            #TODO
            path = bridge_creation(self,times)
            return path
    
    def drift(self, x, t):
        return self.theta*(self.mu-x)
    def diffusion( self, x,t ):
        return self.sigma
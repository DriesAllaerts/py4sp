#!/usr/bin/env python

'''
Module with some new python functions

Author: Dries Allaerts
Date: May 10, 2017
'''

import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.ndimage.filters import uniform_filter1d
from matplotlib.colors import Normalize

def normalised_correlate(a,v,mode):
    a = (a - np.mean(a))/(np.std(a) * len(a))
    v = (v - np.mean(v))/np.std(v)
    return np.correlate(a,v,mode)

def autocorrelation(x):
    mean = x.mean()
    var = np.var(x)
    xp = x-mean
    corr = np.correlate(xp,xp,'full')[len(x)-1:]/var/len(x)
    lags = np.arange(len(corr))
    return corr,lags

def heaviside(x):
    '''
    Heaviside function:
    -------------------
    1.0 for x > 0
    0.5 for x = 0
    0.0 for x < 0
    a tolerance of 1e-12 is allowed
    '''
    return 0.5*(np.sign(np.around(x,12))+1)
def step(x):
    '''
    Step function:
    --------------
    1.0 for x > 0
    0.0 for x <= 0
    a tolerance of 1e-12 is allowed
    '''
    return 1.0 * (np.around(x,12) > 0.0)
def pulse(x,delta):
    '''
    Pulse function:
    ---------------
    1.0 for abs(x) < delta
    0.5 for abs(x) = delta
    0.0 for abs(x) > delta
    '''
    return step(delta-np.abs(x))
def smoothstep(x):
    '''
    Smooth step function:
    ---------------------
    0.0                    for x <= 0
    1/[1+exp(1/(x-1)+1/x)] for 0 < x < 1
    1.0                    for x >= 1
    '''
    with np.errstate(divide='ignore',over='ignore'):
        S = 1.0/(1.0+np.exp(1/(x-1)+1/x)) * ( np.all([x > 0.0,x < 1.0],axis=0) )
    S[x >= 1.0] = 1.0
    return S

def smooth(x, N,mode='nearest'):
    if N>0:
        return uniform_filter1d(x, size=N,mode=mode)
    else:
        return x

def trapzst(u,zst,i1,i2):
    '''Integrate a given profile from i1 to i2,
    with zst the location of the grid cell faces'''
    dz = zst[1:]-zst[:-1]
    integral = np.sum(u[...,i1:i2+1]*dz[i1:i2+1],axis=-1)
    return integral

def plane2image(plane):
    '''Transpose and flip 2D plane to get data in image format'''
    return np.flipud(np.transpose(plane))

def disk_areas(zst,zc,r):
    '''
    For a circle with given height and radius, find the area of the circle
    segment at every grid cell of a given vertical grid

    Parameters
    ----------
    zst: numpy 1D array of Nz+1
        faces of vertical grid cells (=staggered grid)
    zc,r: float
        center and radius of the circle

    Returns (integral)
    ------------------
    areas: numpy 1D array
        areas of circle sigments
    '''
    Nz = zst.shape[0]-1
    areas = np.zeros((Nz))
    for k in range(Nz):
        zl = zst[k]
        zh = zst[k+1]
        if (zl<=zc-r) and np.abs(zh-zc)<r:
            theta = 2*np.arccos(np.abs(zc-zh)/r)
            areas[k] = 0.5*r**2*(theta-np.sin(theta))
        elif np.abs(zh-zc)<r and np.abs(zl-zc)<r:
            if zh<=zc or zl>=zc:
                theta1 = 2*np.arccos(np.abs(zc-zh)/r)
                area1  = 0.5*r**2*(theta1-np.sin(theta1))
                theta2 = 2*np.arccos(np.abs(zc-zl)/r)
                area2  = 0.5*r**2*(theta2-np.sin(theta2))
                areas[k] = np.abs(area1-area2)
            elif zl<zc and zh>zc:
                theta1 = 2*np.arccos(np.abs(zc-zh)/r)
                area1  = 0.5*r**2*(theta1-np.sin(theta1))
                theta2 = 2*np.arccos(np.abs(zc-zl)/r) 
                area2  = 0.5*r**2*(theta2-np.sin(theta2))
                areas[k] = np.pi*r**2-area1-area2
            else:
                print('Error, you should not end up here')
                return 1
        elif (zh>=zc+r) and np.abs(zl-zc)<r:
            theta = 2*np.arccos(np.abs(zc-zl)/r)
            areas[k] = 0.5*r**2*(theta-np.sin(theta))
        #else: point lies outside the circle so do not
    return areas

class progressBar(object):
    def __init__(self,text,N):
        print('    '+text)
        print('    1----------------------------------------100')
        print('     ',end='')

        self.__N = N
        self.__counter = 0

    @property
    def N(self):
        return self.__N

    def incr(self,n):
        dn = round(n/self.N*40)
        ni = dn-self.__counter
        for i in range(ni):
            print('|',end='',flush=True)
        
        self.__counter = dn
        #Reached the end
        if n==self.N:
            print('\n')

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

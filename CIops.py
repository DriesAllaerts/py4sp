#!/usr/bin/env python

'''
Functions to estimate capping inversion parameters

Author: Dries Allaerts
Date: May 22, 2017
'''

import numpy as np
from scipy.optimize import curve_fit

def RZfit(z,th,p0=None):
    '''Estimate parameters of a smooth analytical curve to fit a given
    vertical potential temperature profile

    Parameters
    ----------
    z, th: numpy 1D array
        height and potential temperature profile
    p0: list of length 5
        Initial guess for [a,b,thm,l,dh]

    Returns (CIdata)
    ----------------
    CIdata: dict
        capping inversion parameters
    '''
    popt, pcov = curve_fit(RZmodel,z,th,p0)
    CIdata = {}
    CIdata['a']   = popt[0]
    CIdata['b']   = popt[1]
    CIdata['thm'] = popt[2]
    CIdata['l']   = popt[3]
    CIdata['dh']  = popt[4]
    CIdata['std'] = np.sqrt(np.diag(pcov))

    CIdata['ksi']   = 1.5
    CIdata['c']     = 1.0/(2*CIdata['ksi'])
    CIdata['gamma'] = CIdata['b']/(CIdata['c']*CIdata['dh'])
    CIdata['h0']    = CIdata['l'] - CIdata['dh']/2
    CIdata['h1']    = CIdata['l']
    CIdata['h2']    = CIdata['l'] + CIdata['dh']/2
    CIdata['th00']  = CIdata['a'] + CIdata['thm'] - CIdata['gamma']*CIdata['l']
    CIdata['dth']   = CIdata['a'] + CIdata['b']
    CIdata['dthp']  = CIdata['a']
    CIdata['thfit'] = RZmodel(z,CIdata['a'],CIdata['b'],CIdata['thm'],CIdata['l'],CIdata['dh'])
    return CIdata

def RZmodel(z,a,b,thm,l,dh):
    '''
    Smooth curve representing the vertical potential temperature profile
    of a neutral atmospheric boundary layer with a capping inversion
    Rampanelli & Zardi (2004)

    Parameters
    ----------
    z: numpy 1D array
        height
    a,b,thm,l,dh: float
        fitting parameters

    Returns (th)
    ------------
    th: numpy 1D array
        vertical potential temperature profile
    '''
    ksi = 1.5
    c = 1.0/(2*ksi)
    eta = (z-l)/(c*dh)

    f = (np.tanh(eta)+1.0)/2.0
    with np.errstate(over='ignore'):
        g = (np.log(2*np.cosh(eta))+eta)/2.0
    for i in np.where(np.isinf(g)):
        g[i] = (np.abs(eta)+eta)/2.0

    th = thm + a*f + b*g
    return th


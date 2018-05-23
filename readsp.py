#!/usr/bin/env python

'''
Module for reading setup files

Author: Dries Allaerts
Date: May 10, 2017
'''

import numpy as np
import os

def zmesh(filename, Lref=1.0):
    '''
    read zmesh file

    Parameters
    ----------
    filename: str
        name of the mesh file
    Lref (optional): float
        reference length
        default: 1.0 m

    Returns (Nz,zcc,zst)
    --------------------
    Nz: int
        number of grid cells in the vertical direction
    zcc, zst: 1d numpy array
        cell-centered and staggered vertical grid
    '''
    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist, ZMESH file not found')
        return 1
    
    dummy = np.loadtxt(filename)
    Nz = int((dummy[0]-1)/2)
    zst = dummy[1::2]*Lref
    zcc = dummy[2::2]*Lref
    return Nz, zcc, zst

def ns_setup(filename='NS.setup',Lref=1.0,Uref=1.0):
    '''
    read NS setup file

    Parameters
    ----------
    filename (optional): str
        name of the NS setup file
        default: NS.setup (in current directory)
    Lref, Uref (optional): float
        reference length and velocity
        default: 1.0 m, 1.0 m/s

    Returns (ns_stp)
    ----------------
    ns_stp: dict
        info read from NS.setup stored in dictionary
        keys > case: str
               Nx2, Ny: int
               Lx, Ly: float
               force_mode: int
               force_param: float
               tstart, tstop, dt1, dt2: float
               zmeshf: str
               z0: float
               model_type: int
    '''
    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist, NS setup file not found')
        return 1

    tref = Lref/Uref
    input = []
    with open(filename,'r') as file:
        for line in file:
            input.append(line.rstrip('\r\n').split(' '))

    ns_stp = {}
    ns_stp['case']        = input[1][0]
    ns_stp['Nx2']         = int(input[3][0])
    ns_stp['Ny']          = int(input[3][1])
    ns_stp['Lx']          = float(input[5][0].replace('d','e'))*Lref
    ns_stp['Ly']          = float(input[5][1].replace('d','e'))*Lref
    ns_stp['force_mode']  = int(input[9][0])
    if ns_stp['force_mode']==1:
        scale = 1.0
    elif ns_stp['force_mode']==2:
        scale = Uref**2/Lref
    elif ns_stp['force_mode']==3:
        scale = Uref

    ns_stp['force_param'] = float(input[9][1].replace('d','e'))*scale
    ns_stp['tstart']      = float(input[11][0].replace('d','e'))*tref
    ns_stp['tstop']       = float(input[13][0].replace('d','e'))*tref
    ns_stp['dt1']         = float(input[15][0].replace('d','e'))*tref
    ns_stp['dt2']         = float(input[17][0].replace('d','e'))*tref
    ns_stp['zmeshf']      = input[21][0]
    ns_stp['z0']          = float(input[27][2].replace('d','e'))*Lref
    ns_stp['model_type']  = int(input[29][0])

    return ns_stp

def en_setup(filename='therm.setup',Lref=1.0,Uref=1.0):
    '''
    read EN setup file

    Parameters
    ----------
    filename (optional): str
        name of the EN setup file
        default: therm.setup (in current directory)
    Lref, Uref (optional): float
        reference length and velocity
        default: 1.0 m, 1.0 m/s

    Returns (en_stp)
    ----------------
    en_stp: dict
        info read from therm.setup stored in dictionary
        keys > gravity: float
               thermo: bool
    '''
 
    en_stp = {}
    if not os.path.exists(filename):
        en_st['thermo'] = False
    else:
        input = []
        with open(filename,'r') as file:
            for line in file:
                input.append(line.rstrip('\r\n').split(' '))
    
        en_stp['thermo']    = input[1][0] in ['.TRUE.','.True.','.true.']
        en_stp['gravity']   = float(input[3][0].replace('d','e'))*Uref**2/Lref
    return en_stp

def ek_setup(filename='ekman.setup',Lref=1.0,Uref=1.0):
    '''
    read EK setup file

    Parameters
    ----------
    filename (optional): str
        name of the EK setup file
        default: ekman.setup (in current directory)
    Lref, Uref (optional): float
        reference length and velocity
        default: 1.0 m, 1.0 m/s

    Returns (ek_stp)
    ----------------
    ek_stp: dict
        info read from ekman.setup stored in dictionary
        keys > fc: float
               ekman: bool
    '''
    ek_stp = {}
    if not os.path.exists(filename):
        ek_stp['ekman'] = False
    else:
        input = []
        with open(filename,'r') as file:
            for line in file:
                input.append(line.rstrip('\r\n').split(' '))
        
        ek_stp['ekman'] = input[1][0] in ['.TRUE.','.True.','.true.']
        Ro = float(input[3][0].replace('d','e'))
        ek_stp['fc']   = 1.0/Ro*Uref/Lref
    return ek_stp

def recl_setup(filename='recycling.setup',Lref=1.0,Uref=1.0):
    '''
    read recycling setup file

    Parameters
    ----------
    filename (optional): str
        name of the recycling setup file
        default: recycling.setup (in current directory)
    Lref, Uref (optional): float
        reference length and velocity
        default: 1.0 m, 1.0 m/s

    Returns (recl_stp)
    ----------------
    recl_stp: dict
        info read from recycling.setup stored in dictionary
        keys > non_periodic: bool
               inflow_type: str
               Nx2: int
    '''

    recl_stp = {}
    if not os.path.exists(filename):
        recl_stp['non_periodic'] = False
        recl_stp['inflow_type']  = ''
        return recl_stp
    
    input = []
    with open(filename,'r') as file:
        for line in file:
            input.append(line.rstrip('\r\n').split(' '))

    recl_stp['non_periodic'] = input[1][0] in ['.TRUE.','.True.','.true.']
    recl_stp['inflow_type']  = input[3][0]
    recl_stp['Nx2']          = int(input[5][0])
    recl_stp['Lxfringe']     = (float(input[13][0].replace('d','e')) - float(input[11][0].replace('d','e')))/100.0
    recl_stp['lambda_max']   = float(input[15][0].replace('d','e'))*Uref/Lref
    recl_stp['drise']        = float(input[17][0].replace('d','e'))/100.0
    recl_stp['dfall']        = float(input[19][0].replace('d','e'))/100.0
    return recl_stp

def windfarm_setup(filename='windfarm.setup',Lref=1.0,Uref=1.0):
    '''
    read windfarm setup file

    Parameters
    ----------
    filename (optional): str
        name of the windfarm setup file
        default: windfarm.setup (in current directory)
    Lref, Uref (optional): float
        reference length and velocity
        default: 1.0 m, 1.0 m/s

    Returns (Nt,Ct,windfarm)
    ------------------------
    Nt: int
        number of turbines
    Ct: float
        thrust coefficient
    windfarm: numpy array
        turbine data
    '''
    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist, windfarm setup file not found')
        return 1
   
    tref = Lref/Uref
    with open(filename,'r') as file:
        Nt = int(file.readline().rstrip('\r\n'))
        Ct = float(file.readline().rstrip('\r\n').split()[0])
        data = [] #Empty list
        for line in file:
            data.append([float(i) for i in line.rstrip('\r\n').split()])
        windfarm = np.array(data)
        windfarm[:,0:5] *= Lref #First 5 floats represent a length
        windfarm[:,5] *= 1.0/tref
    return Nt, Ct, windfarm

def PBSout(filename,verbose=True):
    walltime = []
    timestep = []
    total_walltime = None
    with open(filename,'r') as file:
        for line in file:
#            if 'PBS: job killed' in line:
#                break
            if 'RK_time,' in line.rstrip('\r\n').split():
                walltime.append(float(line.rstrip('\r\n').split()[-1]))
                timestep.append(float(line.rstrip('\r\n').split()[-2]))
#            if 'Resource List' in line.rstrip('\r\n').split(':'):
#                node_string = line.rstrip('\r\n').split(',')[1]
#                nodes = int(node_string.replace(':','=').split('=')[1])
#                cores = int(node_string.replace(':','=').split('=')[-1])
            if 'nodes' in  line.rstrip('\r\n').split(': '):
                nodes = int(line.rstrip('\r\n').split(': ')[-1])
            if 'procs' in  line.rstrip('\r\n').split(': '):
                cores = int(int(line.rstrip('\r\n').split(': ')[-1])/nodes)
            if 'Resources Used' in line.rstrip('\r\n').split(':'):
                walltime_string = line.rstrip('\r\n').split(',')[-1]
                hours = float(walltime_string.replace(':','=').split('=')[1])
                minutes = float(walltime_string.replace(':','=').split('=')[2])
                seconds = float(walltime_string.replace(':','=').split('=')[3])
                total_walltime = hours+minutes/60.+seconds/3600.
    if not total_walltime:
        total_walltime=np.sum(walltime)/3600.
    N = len(walltime)
    if verbose:
        print('Average time step for',N,'samples is',np.mean(np.array(timestep)))
        print('Average walltime per time step for',N,'samples is',np.mean(np.array(walltime)))
        print('Number of cores:',nodes,'x',cores)
        print('Total wall time:',total_walltime,'h')
    return nodes,cores,total_walltime,N,np.array(walltime)


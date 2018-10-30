#!/usr/bin/env python

'''
Module for loading data

Author: Dries Allaerts
Date: May 10, 2017
'''

import numpy as np
import os
from py4sp import fieldsp

def NSpost(filename='NS_post1.dat',case='BL',model_type=1,Lref=1.0,Uref=1.0,unique_time=True):
    '''
    Load NS_post1.dat file
    '''
    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    #Scaling    
    tref = Lref/Uref

    #Read file
    dummy = np.loadtxt(filename,skiprows=1)
    NSp = {}
    if unique_time:
        NSp['t'],indices = np.unique(dummy[:,0]*tref,return_index=True)
    else:
        NSp['t'] = dummy[:,0]*tref
        indices = range(dummy.shape[0])

    NSp['U'] = dummy[indices,1]*Uref
    NSp['V'] = dummy[indices,2]*Uref
    NSp['W'] = dummy[indices,3]*Uref
    NSp['Etot'] = dummy[indices,4]*Uref**2
    NSp['Eturb'] = dummy[indices,5]*Uref**2
    if case=='channel':
        if model_type==0:
            NSp['sf_min'] = dummy[indices,6]*Uref**2
            NSp['sf_plus'] = dummy[indices,7]*Uref**2
            NSp['Re'] = dummy[indices,8]
            NSp['utau'] = dummy[indices,9]*Uref
        else:
            NSp['sf_lmin'] = dummy[indices,6]*Uref**2
            NSp['sf_lplus'] = dummy[indices,7]*Uref**2
            NSp['sf_tmin'] = dummy[indices,8]*Uref**2
            NSp['sf_tplus'] = dummy[indices,9]*Uref**2
            NSp['Re'] = dummy[indices,10]
            NSp['utau'] = dummy[indices,11]*Uref
    else:
        if model_type==0:
            NSp['sf_mean_x'] = dummy[indices,6]*Uref**2 # based on integral of right-hand side stress terms
            NSp['sf_mean_y'] = dummy[indices,7]*Uref**2 # based on integral of right-hand side stress terms
            NSp['sf_Um'] = dummy[indices,8]*Uref**2 # based on log law an <U>[0]
        elif model_type==1:
            NSp['sf_mean_x'] = dummy[indices,6]*Uref**2 # based on integral of right-hand side stress terms
            NSp['sf_mean_y'] = dummy[indices,7]*Uref**2 # based on integral of right-hand side stress terms
            NSp['sf_Um'] = dummy[indices,8]*Uref**2 # based on log law an <U>[0]
            NSp['sf_sgs'] = dummy[indices,9]*Uref**2
        elif model_type==-2:
            NSp['Esgs'] = dummy[indices,6]*Uref**2
            NSp['sf_Um'] = dummy[indices,7]*Uref**2
        else:
            print('Error: model_type unknown')
            return 1

    return NSp

def ENpost(filename='en_post1.dat',Lref=1.0,Uref=1.0,Tref=288.15,unique_time=True):
    '''Load en_post1.dat file
    Current implementation not for channel cases'''
    
    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    #Scaling
    tref = Lref/Uref

    dummy = np.loadtxt(filename,skiprows=1)
    ENp = {}
    if unique_time:
        ENp['t'],indices = np.unique(dummy[:,0]*tref,return_index=True)
    else:
        ENp['t'] = dummy[:,0]*tref
        indices = range(dummy.shape[0])

    ENp['th'] = dummy[indices,1]*Tref
    ENp['th2'] = dummy[indices,2]*Tref**2
    ENp['ttau'] = dummy[indices,3]*Tref
    ENp['dzeta'] = dummy[indices,4]
    ENp['hBL'] = dummy[indices,5]*Lref
    ENp['T2'] = dummy[indices,6]*Tref
    return ENp

def EKpost(filename='ek_post.dat',Lref=1.0,Uref=1.0,unique_time=True):
    '''Load ek_post.dat file
    Current implementation not for channel cases'''
    
    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    #Scaling
    tref = Lref/Uref

    dummy = np.loadtxt(filename,skiprows=1)
    EKp = {}
    if unique_time:
        EKp['t'],indices = np.unique(dummy[:,0]*tref,return_index=True)
    else:
        EKp['t'] = dummy[:,0]*tref
        indices = range(dummy.shape[0])

    EKp['alpha']     = dummy[indices,1]
    EKp['phi']       = dummy[indices,2]
    EKp['omega_m']   = dummy[indices,3]/tref
    EKp['omega_eff'] = dummy[indices,4]/tref
    return EKp

def BLfield(filename='BL_field.dat',Nl=3,Uref=1.0,Tref=1.0,grid=None):
    '''Load BLfield file (binary)'''

    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    BL = {}
    with open(filename,'rb') as binfile:
        BL['time']  = np.asscalar(np.fromfile(binfile, dtype=np.float64, count=1))
        BL['Lx']    = np.asscalar(np.fromfile(binfile, dtype=np.float64, count=1))
        BL['Ly']    = np.asscalar(np.fromfile(binfile, dtype=np.float64, count=1))
        BL['Nx2']   = np.asscalar(np.fromfile(binfile, dtype=np.int32, count=1))
        BL['Ny']    = np.asscalar(np.fromfile(binfile, dtype=np.int32, count=1))
        BL['Nz']    = np.asscalar(np.fromfile(binfile, dtype=np.int32, count=1))
        BL['alpha'] = np.asscalar(np.fromfile(binfile, dtype=np.float64, count=1))
        if Nl>3:
            BL['thetaground'] = np.asscalar(np.fromfile(binfile, dtype=np.float64, count=1))
        #Determine size
        BL['Nx'] = int(BL['Nx2']/2+1) #Plus 1 for defunct mode
        amount  = BL['Nx']*BL['Ny']*BL['Nz']
        amount2 = BL['Nx']*BL['Ny']*(BL['Nz']-1)
        shape   = (BL['Nx'],BL['Ny'],BL['Nz'])
        shape2  = (BL['Nx'],BL['Ny'],BL['Nz']-1)
        #Don't use a dummy array for memory efficiency
        BL['u'] = fieldsp.cfield(np.fromfile(binfile,dtype=np.complex128,count=amount).reshape(shape,order='F')*Uref,grid)
        BL['v'] = fieldsp.cfield(np.fromfile(binfile,dtype=np.complex128,count=amount).reshape(shape,order='F')*Uref,grid)
        BL['w'] = fieldsp.cfield(np.fromfile(binfile,dtype=np.complex128,count=amount2).reshape(shape2,order='F')*Uref,grid)
        if Nl>3:
            BL['th'] = fieldsp.cfield(np.fromfile(binfile,dtype=np.complex128,count=amount).reshape(shape,order='F')*Tref,grid)
        if Nl>4:
            BL['tke'] = fieldsp.cfield(np.fromfile(binfile,dtype=np.complex128,count=amount).reshape(shape,order='F')*Uref**2,grid)
        BL['kx'] = np.array([(i)/BL['Lx']*(2*np.pi) for i in range(BL['Nx'])])
        BL['ky'] = np.array([(i)/BL['Ly']*(2*np.pi) for i in range(int(-BL['Ny']/2+1),int(BL['Ny']/2+1))])
    return BL

def BLfieldstat(filename='BL_fieldstat.dat',Nl=11,Uref=1.0,**kwargs):
    '''Load BLfieldstat file (single binary)'''

    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    if 'grid' in kwargs:
        grid = kwargs['grid']
        Nx2 = grid.Nx2
        Ny  = grid.Ny
        Nz  = grid.Nz
    elif all([i in kwargs for i in ['Nx2','Ny','Nz']]):
        grid = None
        Nx2 = kwargs['Nx2']
        Ny  = kwargs['Ny']
        Nz  = kwargs['Nz']
    else:
        print('Error: Cannot determine the dimensions of BLfieldstat')
        return 1

    stat = {}
    with open(filename,'rb') as binfile:
        stat['nsamp'] = np.fromfile(binfile, dtype=np.int32, count=1)
        stat['time_interv'] = np.fromfile(binfile, dtype=np.float32, count=1)
        stat['time_incurr'] = np.fromfile(binfile, dtype=np.float32, count=1)
        dummy = np.fromfile(binfile, dtype=np.float64, count=Nx2*Ny*Nz*Nl)
    shape = (Nx2,Ny,Nz,Nl)
    dummy = dummy.reshape(shape, order='F')
    stat['u']   = fieldsp.rfield(dummy[:,:,:,0]*Uref,grid)
    stat['v']   = fieldsp.rfield(dummy[:,:,:,1]*Uref,grid)
    stat['w']   = fieldsp.rfield(dummy[:,:,:,2]*Uref,grid)
    stat['uu']  = fieldsp.rfield(dummy[:,:,:,3]*Uref**2,grid)
    stat['vv']  = fieldsp.rfield(dummy[:,:,:,4]*Uref**2,grid)
    stat['ww']  = fieldsp.rfield(dummy[:,:,:,5]*Uref**2,grid)
    stat['uv']  = fieldsp.rfield(dummy[:,:,:,6]*Uref**2,grid)
    stat['uw']  = fieldsp.rfield(dummy[:,:,:,7]*Uref**2,grid)
    stat['vw']  = fieldsp.rfield(dummy[:,:,:,8]*Uref**2,grid)
    stat['p']   = fieldsp.rfield(dummy[:,:,:,9]*Uref**2,grid)
    stat['wst'] = fieldsp.rfield(dummy[:,:,:,10]*Uref,grid)

    #names = ['u', 'v', 'w', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw', 'p', 'wst']
    #for i in range(Nl):
    #    stat[names[i]] = fieldsp.rfield(dummy[:,:,:,i],grid)
    return stat

def BLsurfacestat(filename='BL_surface_stat.dat',Uref=1.0,Tref=1.0,**kwargs):
    '''Load BL_surface_stat file'''

    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    if 'grid' in kwargs:
        grid = kwargs['grid']
        Nx2 = grid.Nx2
        Ny  = grid.Ny
    elif all([i in kwargs for i in ['Nx2','Ny']]):
        grid = None
        Nx2 = kwargs['Nx2']
        Ny  = kwargs['Ny']
    else:
        print('Error: Cannot determine the dimensions of BLsurfacestat')
        return 1

    stat = {}
    Nl = 3
    with open(filename,'rb') as binfile:
        stat['nsamp'] = np.fromfile(binfile, dtype=np.int32, count=1)
        dummy = np.fromfile(binfile, dtype=np.float64, count=Nx2*Ny*Nl)
    shape = (Nx2,Ny,Nl)
    dummy = dummy.reshape(shape, order='F')
    stat['tw1']  = fieldsp.rfield(dummy[:,:,0]*Uref**2,grid)
    stat['tw2']  = fieldsp.rfield(dummy[:,:,1]*Uref**2,grid)
    stat['qw']   = fieldsp.rfield(dummy[:,:,2]*Uref*Tref,grid)

    return stat

def BLfieldstat_distributed(filename='BL_fieldstat_01.dat',scale=1.0,**kwargs):
    '''Load a particular field from a distributed BLfieldstat file (binary)'''

    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    if 'grid' in kwargs:
        grid = kwargs['grid']
        Nx2 = grid.Nx2
        Ny  = grid.Ny
        Nz  = grid.Nz
    elif all([i in kwargs for i in ['Nx2','Ny','Nz']]):
        grid = None
        Nx2 = kwargs['Nx2']
        Ny  = kwargs['Ny']
        Nz  = kwargs['Nz']
    else:
        print('Error: Cannot determine the dimensions of BLfieldstat')
        return 1

    metadata = {}
    with open(filename,'rb') as binfile:
        metadata['nsamp'] = np.fromfile(binfile, dtype=np.int32, count=1)
        metadata['time_interv'] = np.fromfile(binfile, dtype=np.float32, count=1)
        metadata['time_incurr'] = np.fromfile(binfile, dtype=np.float32, count=1)
        #dummy = np.fromfile(binfile, dtype=np.float64, count=Nx2*Ny*Nz)
        dummy = np.fromfile(binfile, dtype=np.float64, count=-1)
    shape = (Nx2,Ny,Nz)
    dummy = dummy.reshape(shape, order='F')
    return fieldsp.rfield(dummy*scale,grid), metadata

def BLtstatcc(filename='BL_tstatcc.dat',Lref=1.0,Uref=1.0):
    '''Load BL_tstatcc file'''
 
    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    with open(filename,'r') as file:
        input = file.readline().rstrip('\r\n').split()

    data = {}
    data['Nsample'] = int(input[4])
    data['utau'] = float(input[9])*Uref
    #dummy = np.loadtxt(filename,skiprows=2) #Has some problems when data is corrupted (e.g. some files contain ...-270 instead of ....E-270 due to a bug in SP-Wind)
    dummy = np.genfromtxt(filename,skip_header=2,loose=True)
    data['z']  = dummy[:,0]*Lref
    data['u']  = dummy[:,1]*Uref
    data['v']  = dummy[:,2]*Uref
    data['w']  = dummy[:,3]*Uref
    data['uu'] = dummy[:,4]*Uref**2
    data['vv'] = dummy[:,5]*Uref**2
    data['ww'] = dummy[:,6]*Uref**2
    data['uv'] = dummy[:,7]*Uref**2
    data['uw'] = dummy[:,8]*Uref**2
    data['vw'] = dummy[:,9]*Uref**2
    return data

def BLtstatst(filename='BL_tstatst.dat',Lref=1.0,Uref=1.0):
    '''Load BL_tstatst file'''
 
    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    with open(filename,'r') as file:
        input = file.readline().rstrip('\r\n').split()

    data = {}
    data['Nsample'] = int(input[4])
    data['utau'] = float(input[9])*Uref
    #dummy = np.loadtxt(filename,skiprows=2) #Has some problems when data is corrupted
    dummy = np.genfromtxt(filename,skip_header=2,loose=True)
    data['zst']   = dummy[:,0]*Lref
    data['dudz']  = dummy[:,1]*Uref/Lref
    data['dvdz']  = dummy[:,2]*Uref/Lref
    data['s13']   = dummy[:,3]*Uref**2
    data['s23']  = dummy[:,4]*Uref**2
    return data

def BLtstatall(ccfilename='BL_tstatcc.dat',stfilename='BL_tstatst.dat',Lref=1.0,Uref=1.0):
    '''Load BL_tstatcc and BL_tstatst files and add some additional data'''
    datacc = BLtstatcc(ccfilename,Lref,Uref)
    datast = BLtstatst(stfilename,Lref,Uref)
    data = dict(datacc,**datast)

    return data

def ENtstatcc(filename='en_tstatcc.dat',Lref=1.0,Uref=1.0,Tref=1.0):
    '''Load en_tstatcc file'''
 
    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1

    with open(filename,'r') as file:
        input = file.readline().rstrip('\r\n').split()

    data = {}
    data['Nsample'] = int(input[4])
    dummy = np.loadtxt(filename,skiprows=2)
    data['z']  = dummy[:,0]*Lref
    data['w']  = dummy[:,1]*Uref
    data['th']  = dummy[:,2]*Tref
    data['th2']  = dummy[:,3]*Tref**2
    data['wth'] = dummy[:,4]*Uref*Tref
    return data

def Windpower(filename='Windpower.dat',Lref=1.0,Uref=1.0,unique_time=True):
    '''Load Windpower file'''

    if not os.path.exists(filename):
        print('Error: '+filename+' does not exist')
        return 1
    
    tref = Lref/Uref
    dummy = np.loadtxt(filename,comments='%')
    if unique_time:
        time,indices = np.unique(dummy[:,0]*tref,return_index=True)
    else:
        time = dummy[:,0]*tref
        indices = range(dummy.shape[0])

    force = -dummy[indices,1::2]*Uref**2*Lref**2
    power = -dummy[indices,2::2]*Uref**3*Lref**2
    return time, force, power

def ztfile(filename,Nz,Nt,scale=1.0):
    '''Load a zt file (binary)'''
    with open(filename,'rb') as binfile:
        dummy = np.fromfile(binfile, dtype=np.float64, count=Nz*Nt)
        shape = (Nz,Nt)
        field = dummy.reshape(shape, order='F')*scale
    return field

def tfile(filename,tref=1.0):
    '''Load time_array file (ascii)'''
    return tref*np.loadtxt(filename,skiprows=1)

def tfile_bin(filename,tref=1.0):
    '''Load time_array file (binary)'''
    with open(filename,'rb') as binfile:
        Nt = np.asscalar(np.fromfile(binfile, dtype=np.int32, count=1))
        ts = np.fromfile(binfile, dtype=np.float64, count=Nt)*tref
    return ts

def CIfile_matlab(filename):
    '''Load file with Capping Inversion data
    (old data format generated with matlab)'''
    dummy = np.loadtxt(filename,delimiter=',')
    CI = {}
    CI['thm']   = dummy[:,1]
    CI['dh']    = dummy[:,2]
    CI['dth']   = dummy[:,3]
    CI['h0']    = dummy[:,4]
    CI['h1']    = dummy[:,5]
    CI['h2']    = dummy[:,6]
    CI['th00']  = dummy[:,7]
    CI['gamma'] = dummy[:,8]
    CI['sse']   = dummy[:,9]
    CI['rmse']  = dummy[:,10]
    CI['a']     = dummy[:,11]
    CI['b']     = dummy[:,12]
    CI['c']     = dummy[:,13]
    CI['l']     = dummy[:,14]
    return CI

def point_measurement(filename,Uref=1.0,Lref=1.0):
    dummy = np.loadtxt(filename)
    data = {}
    data['t'] = dummy[:,0]*Lref/Uref
    data['u'] = dummy[:,1]*Uref
    data['v'] = dummy[:,2]*Uref
    data['w'] = dummy[:,3]*Uref
    return data

def image(filename):
    with open(filename,'rb') as binfile:
        m = np.asscalar(np.fromfile(binfile, dtype=np.int32, count=1))
        n = np.asscalar(np.fromfile(binfile, dtype=np.int32, count=1))
        data = np.fromfile(binfile,dtype=np.float64,count=m*n).reshape((m,n))
    return data

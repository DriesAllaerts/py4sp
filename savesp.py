#!/usr/bin/env python

'''
Module for saving data to disk

Author: Dries Allaerts
Date: May 17, 2017
'''

import numpy as np
import os
from py4sp import fieldsp

def BLfield(BL,filename='BL_field.dat',Nl=3):
    '''Save BLfield file (binary)'''

    with open(filename,'wb') as binfile:
        np.array([BL['time']]).astype(dtype=np.float64).tofile(binfile)
        np.array([BL['Lx']]).astype(dtype=np.float64).tofile(binfile)
        np.array([BL['Ly']]).astype(dtype=np.float64).tofile(binfile)
        np.array([BL['Nx2']]).astype(dtype=np.int32).tofile(binfile)
        np.array([BL['Ny']]).astype(dtype=np.int32).tofile(binfile)
        np.array([BL['Nz']]).astype(dtype=np.int32).tofile(binfile)
        np.array([BL['alpha']]).astype(dtype=np.float64).tofile(binfile)
        if Nl>3:
            np.array([BL['thetaground']]).astype(dtype=np.float64).tofile(binfile)
        BL['u'].value.reshape((BL['u'].size,1), order='F').tofile(binfile)
        BL['v'].value.reshape((BL['v'].size,1), order='F').tofile(binfile)
        BL['w'].value.reshape((BL['w'].size,1), order='F').tofile(binfile)
        if Nl>3:
            BL['th'].value.reshape((BL['th'].size,1), order='F').tofile(binfile)
        if Nl>4:
            BL['tke'].value.reshape((BL['tke'].size,1), order='F').tofile(binfile)

def image(filename,data):
    with open(filename,'wb') as binfile:
        np.array([data.shape[0]]).astype(dtype=np.int32).tofile(binfile)
        np.array([data.shape[1]]).astype(dtype=np.int32).tofile(binfile)
        data.ravel().tofile(binfile)

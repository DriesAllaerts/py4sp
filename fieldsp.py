#!/usr/bin/env python

'''
Module defining three-dimensional scalar and vector fields

Author: Dries Allaerts
Date: May 11, 2017
'''

import numpy as np
import os
import matplotlib.pyplot as plt
from py4sp import mypy

class cfield(object):
    def __init__(self,value,grid=None):
        self.__value = value
        self.__grid = grid
        self.__shape = value.shape
        self.__size  = value.size

    @property
    def value(self):
        return self.__value

    @property
    def grid(self):
        return self.__grid

    @property
    def shape(self):
        return self.__shape

    @property
    def size(self):
        return self.__size

    def c2r(self):
        res = np.fft.ifft(self.value,axis=1)
        res = np.fft.irfft(res,axis=0)
        Nx2 = (self.shape[0]-1)*2
        Ny  = self.shape[1]
        return rfield(res*Nx2*Ny,self.grid)

class rfield(object):
    def __init__(self,value,grid=None):
        self.__value = value
        self.__grid = grid
        self.__shape = value.shape

    @property
    def value(self):
        return self.__value

    @property
    def grid(self):
        return self.__grid

    @property
    def shape(self):
        return self.__shape

    def r2c(self):
        res = np.fft.rfft(self.value,axis=0)
        res = np.fft.fft(res,axis=1)
        Nx2 = self.shape[0]
        Ny  = self.shape[1]
        return cfield(res/Nx2/Ny,self.grid)
    
    def xyplane(self,z):
        k = np.max(np.where(self.grid.zcc<=z))
        plane = self.value[:,:,k] + (self.value[:,:,k+1]-self.value[:,:,k])/(self.grid.zcc[k+1]-self.grid.zcc[k])*(z-self.grid.zcc[k])
        return plane

    def xzplane(self,y):
        j = np.max(np.where(self.grid.ys<=y))
        plane = self.value[:,j,:] + (self.value[:,j+1,:]-self.value[:,j,:])/(self.grid.ys[j+1]-self.grid.ys[j])*(y-self.grid.ys[j])
        return plane

    def yzplane(self,x):
        i = np.max(np.where(self.grid.xs<=x))
        plane = self.value[i,:,:] + (self.value[i+1,:,:]-self.value[i,:,:])/(self.grid.xs[i+1]-self.grid.xs[i])*(x-self.grid.xs[i])
        return plane

    def topview(self,z):
        #Using pcolormesh
        #plane = self.xyplane(z)
        #plt.figure()
        #plt.pcolormesh(self.grid.xs,self.grid.ys,np.transpose(plane))
        #plt.axis('equal')
        #plt.colorbar()
        #plt.show(block=False)

        #Alternative: use imshow
        image = mypy.plane2image(self.xyplane(z))
        plt.figure()
        plt.imshow(image,extent=[.0,self.grid.Lx,.0,self.grid.Ly])
        plt.colorbar()
        plt.show(block=False)

    def sideview(self,y):
        plane = self.xzplane(y)
        plt.figure()
        plt.pcolormesh(self.grid.xs,self.grid.zcc,np.transpose(plane))
        plt.axis('equal')
        plt.colorbar()
        plt.show(block=False)

    def frontview(self,x):
        plane = self.yzplane(x)
        plt.figure()
        plt.pcolormesh(self.grid.ys,self.grid.zcc,np.transpose(plane))
        plt.axis('equal')
        plt.colorbar()
        plt.show(block=False)

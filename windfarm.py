#!/usr/bin/env python

'''
Module defining windfarm class

Author: Dries Allaerts
Date: May 10, 2017
'''

import numpy as np
import os
import textwrap
import matplotlib.pyplot as plt
from py4sp import readsp
from py4sp import loadsp
from py4sp import fieldsp


class Windfarm(object):
    '''Wind farm data structure'''

    def __init__(self,path=None,Lref=1.0,Uref=1.0,**kwargs):
        if not path:
            path = os.getcwd()

        if 'setupfile' in kwargs:
            setupfilename = kwargs['setupfile']
        else:
            setupfilename = 'windfarm.setup'

        if 'datafile' in kwargs:
            datafilename = kwargs['datafile']
        else:
            datafilename = 'Windpower.dat'

        #Read setup file
        setupfile = os.path.join(path,setupfilename)
        self.__Nt, self.__Ctp, turb_data = readsp.windfarm_setup(setupfile,Lref=Lref,Uref=Uref)
        self.__x = turb_data[:,0]
        self.__y = turb_data[:,1]
        self.__z = turb_data[:,2]
        self.__r = turb_data[:,4]

        #Check whether windfarm is structured (current only works well for aligned)
        self.__Nxt = np.unique(self.x).shape[0]
        self.__Nyt = np.unique(self.y).shape[0]
        if not self.Nxt*self.Nyt==self.Nt:
            self.__Nxt = None
            self.__Nyt = None
            self.__sx  = 0.0
            self.__sy  = 0.0
        else:
            self.__sx = (self.x[self.Nyt]-self.x[0])/(2*self.r[0])
            self.__sy = (self.y[1]-self.y[0])/(2*self.r[0])

        #Read data file
        datafile = os.path.join(path,datafilename)
        if os.path.exists(datafile):
            self.__t, self.__F, self.__P = loadsp.Windpower(datafile,Lref=Lref,Uref=Uref)
        else:
            self.__t = np.array(())
            self.__F = np.array(())
            self.__P = np.array(())

    @property
    def Nt(self):
        return self.__Nt

    @property
    def Nxt(self):
        return self.__Nxt

    @property
    def Nyt(self):
        return self.__Nyt

    @property
    def sx(self):
        return self.__sx

    @property
    def sy(self):
        return self.__sy

    @property
    def Ctp(self):
        return self.__Ctp

    @property
    def Ct(self):
        return 16.0*self.Ctp/(self.Ctp+4.0)**2
    
    @property
    def x(self):
        return self.__x

    @property
    def y(self):
        return self.__y

    @property
    def z(self):
        return self.__z

    @property
    def r(self):
        return self.__r
    
    @property
    def t(self):
        if self.__t.shape[0]!=0:
            return self.__t
        else:
            print('Error: windfarm data is not loaded')
            return 1

    @property
    def F(self):
        if self.__F.shape[0]!=0:
            return self.__F
        else:
            print('Error: windfarm data is not loaded')
            return 1

    @property
    def P(self):
        if self.__P.shape[0]!=0:
            return self.__P
        else:
            print('Error: windfarm data is not loaded')
            return 1

    def __str__(self):
        windfarm_str='''\
            Wind farm
            ****************************
            Number of turbines: {Nt}
            Thrust coefficient: {Ctp}'''.format(Nt=self.Nt,Ctp=self.Ctp)
        return textwrap.dedent(windfarm_str)

    def Pavg(self,tstart=None,tend=None):
        if tstart:
            istart = np.max(np.where(self.t<=tstart))
        else:
            istart = 0

        if tend:
            iend = np.max(np.where(self.t<=tend))+1
        else:
            iend = len(self.t)

        return np.mean(self.P[istart:iend,:],0)

    def eta_wake(self,tstart=None,tend=None):
        P1   = np.mean(self.Pavg(tstart,tend)[0:self.Nyt])
        Ptot = np.sum(self.Pavg(tstart,tend))
        return Ptot/(self.Nt*P1)

    def Favg(self,tstart=None,tend=None):
        if tstart:
            istart = np.max(np.where(self.t<=tstart))
        else:
            istart = 0

        if tend:
            iend = np.max(np.where(self.t<=tend))+1
        else:
            iend = len(self.t)

        return np.mean(self.F[istart:iend,:],0)

    def plot_layout(self):
        plt.figure()
        for turb in range(self.Nt):
            plt.plot(self.x,self.y,'ok')
        plt.show(block=False)

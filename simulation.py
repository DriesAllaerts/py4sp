#!/usr/bin/env python

'''
Simulation module

Module defining various data structures
'''
__author__ = "Dries Allaerts"
__date__ = "May 10, 2017"

import numpy as np
import os
import textwrap
from py4sp import readsp
from py4sp import loadsp
from py4sp import fieldsp
from py4sp import windfarm as wf

class Simulation(object):
    '''
    data structure containing all data from a particular simulation
    '''
    
    def __init__(self,name=None,path=None,**kwargs):
        '''
        Parameters
        ----------
        name (optional): str
            name of the simulation (e.g. used for labeling and saving)
            default: current directory name
        path (optional): str
            path where setup files and data are read from
            default: current directory
        Lref, Uref, Tref (optional): float
            scaling values to convert SP-Wind data to SI units
            default: 1000 m, 10 m/s, 288.15 K
        verbose (optional): bool
            flag to print information
            default: True
        '''

        if path:
            self.__path = path
        else:
            self.__path = os.getcwd()

        if name:
            self.__name = name
        else:
            self.__name = os.path.basename(self.path)

        #Scaling
        if 'Lref' in kwargs:
            self.__Lref = kwargs['Lref']
        else:
            self.__Lref = 1000.0

        if 'Uref' in kwargs:
            self.__Uref = kwargs['Uref']
        else:
            self.__Uref = 10.0

        if 'Tref' in kwargs:
            self.__Tref = kwargs['Tref']
        else:
            self.__Tref = 288.15

        self.__tref = self.Lref/self.Uref

        #Read NS.setup
        nsfile = os.path.join(self.path,'NS.setup')
        ns_stp = readsp.ns_setup(nsfile,Lref=self.Lref,Uref=self.Uref)
        self.__case = ns_stp['case']
        self.__model_type = ns_stp['model_type']
        self.__grid = Grid(ns_stp['Nx2'],
                         ns_stp['Ny'],
                         ns_stp['Lx'],
                         ns_stp['Ly'],
                         os.path.join(self.path,ns_stp['zmeshf']),
                         Lref=self.Lref)
        
        if self.case=='channel':
            print('Warning: py4sp can not yet handle channel flow cases')

        #Read recycling.setup
        self.__grid_prec = None
        self.__Lxfringe = 0
        self.__lambda_max = 0
        reclfile = os.path.join(self.path,'recycling.setup')
        recl_stp = readsp.recl_setup(reclfile,Lref=self.Lref,Uref=self.Uref)
        if recl_stp['non_periodic'] and recl_stp['inflow_type']=='concurrent_precursor':
            self.__Lxfringe = recl_stp['Lxfringe']*self.grid.Lx
            self.__lambda_max = recl_stp['lambda_max']
            self.__drise = recl_stp['drise']*self.grid.Lx
            self.__dfall = recl_stp['dfall']*self.grid.Lx
            if not (recl_stp['Nx2']==0 or recl_stp['Nx2']==self.grid.Nx2):
                Lx_prec = self.grid.Lx/self.grid.Nx2*recl_stp['Nx2']
                self.__grid_prec = Grid(recl_stp['Nx2'],
                                        ns_stp['Ny'],
                                        Lx_prec,
                                        ns_stp['Ly'],
                                        os.path.join(self.path,ns_stp['zmeshf']),
                                        Lref=self.Lref)
            else:
                self.__grid_prec = self.grid

        #Read windfarm.setup
        self.__wf = None
        if self.case=='Windfarm':
            self.__wf = wf.Windfarm(self.path,Lref=self.Lref,Uref=self.Uref)

        #Generate ABL object
        self.__abl = ABL(self.path,
                Lref=self.Lref,
                Uref=self.Uref,
                Tref=self.Tref)

        #Prepare data structures
        self.__NSp = None
        self.__ENp = None
        self.__fieldstat = None

        #If verbose, print summary of simulation
        if 'verbose' not in kwargs:
            verbose = True
        else:
            verbose = kwargs['verbose']

        if verbose:
            print(self)

    
    @property
    def name(self):
        '''name of the simulation'''
        return self.__name

    @property
    def path(self):
        '''path were setup files and data are read from'''
        return self.__path

    @property
    def case(self):
        '''case type from SP-Wind: channel, BL or Windfarm'''
        return self.__case

    @property
    def Lref(self):
        '''reference length'''
        return self.__Lref

    @property
    def Uref(self):
        '''reference velocity'''
        return self.__Uref

    @property
    def Tref(self):
        '''reference temperature'''
        return self.__Tref

    @property
    def tref(self):
        '''reference time'''
        return self.__tref

    @property
    def model_type(self):
        '''model type: DNS, Smagorinsky, TKE model'''
        return self.__model_type

    @property
    def grid(self):
        '''grid data structure'''
        return self.__grid

    @property
    def grid_prec(self):
        '''
        grid data structure of precursor simulation
        if not defined: None
        '''
        return self.__grid_prec

    @property
    def Lxfringe(self):
        '''
        fringe region length
        if not defined: zero
        '''
        return self.__Lxfringe

    @property
    def lambda_max(self):
        '''
        fringe region damping coefficient
        '''
        return self.__lambda_max

    @property
    def drise(self):
        '''
        drise of fringe region weighting function
        '''
        return self.__drise

    @property
    def dfall(self):
        '''
        dfall of fringe region weighting function
        '''
        return self.__dfall

    @property
    def wf(self):
        '''
        wind farm data structure
        if not defined: None
            '''
        return self.__wf

    @property
    def abl(self):
        '''abl data structure'''
        return self.__abl

    @property
    def NSp(self):
        '''
        dictionary containing data loaded from NS_post1.dat
        if not defined: None
        '''
        return self.__NSp
    
    @property
    def ENp(self):
        '''
        dictionary containing data loaded from en_post1.dat
        if not defined: None
        '''
        return self.__ENp

    @property
    def fieldstat(self):
        '''
        dictionary containing data loaded from BL_fieldstat.dat
        if not defined: None
        '''
        return self.__fieldstat
    
    #Methods for loading data
    def load_NSp(self):
        '''load data from NS_post1.dat into NSp'''
        filename = os.path.join(self.path,'NS_post1.dat')
        self.__NSp = loadsp.NSpost(filename,self.case,self.model_type,Lref=self.Lref,Uref=self.Uref)

    def load_ENp(self):
        '''load data from EN_post1.dat into ENp'''
        filename = os.path.join(self.path,'en_post1.dat')
        self.__ENp = loadsp.ENpost(filename,Lref=self.Lref,Uref=self.Uref,Tref=self.Tref)

    def load_BLfieldstat(self):
        '''
        load data from BL_fieldstat.dat into fieldstat
        (only load when data is not yet loaded to save time)
        '''
        if not self.fieldstat:
            filename = os.path.join(self.path,'BL_fieldstat.dat')
            self.__fieldstat = loadsp.BLfieldstat(filename,Uref=self.Uref,grid=self.grid)
            M = np.sqrt(self.fieldstat['u'].value**2+self.fieldstat['v'].value**2)
            phi = np.arctan(self.fieldstat['v'].value/self.fieldstat['u'].value)*180/np.pi
            self.__fieldstat['M'] = fieldsp.rfield(M,self.grid)
            self.__fieldstat['phi'] = fieldsp.rfield(phi,self.grid)

    def unload_BLfieldstat(self):
        '''unload fieldstat to clean memory'''
        self.__fieldstat = None
    
    def __str__(self):
        sim_str = '''\
            %%%%%%%%%%%%%%%%%%%%%%
            % SP-WIND SIMULATION %
            %%%%%%%%%%%%%%%%%%%%%%
            Name: {name}
            Case: {case}
            Precursor: {prec}
            Path: {path}'''.format(name=self.name,
                                case=self.case,
                                prec=bool(self.grid_prec),
                                path=self.path)
        scaling_str = '''\
            Scaling: Lref = {Lref:.1f} m
                     Uref = {Uref:.2f} m/s
                     tref = {tref:.1f} s
                     Tref = {Tref:.2f} K'''.format(Lref=self.Lref,
        Uref=self.Uref,
        tref=self.tref,
        Tref=self.Tref)
        description = textwrap.dedent(sim_str)+'\n'+textwrap.dedent(scaling_str)
        description += '\n\n'+str(self.grid)
        if self.grid_prec:
            description += '\n\n(Precursor) '+str(self.grid_prec)
            description += '\nFringe region: {Lxfringe:.1f}'.format(Lxfringe=self.Lxfringe)
        description += '\n\n'+str(self.abl)
        if self.wf:
            description += '\n\n'+str(self.wf)
        if any([self.NSp,self.ENp,self.fieldstat]):
            header_str = 'Data\n****************************'
            description += '\n\n'+header_str
            if self.NSp:
                key_str = ', '.join(self.NSp.keys())
                data_str = 'NSpost (Nt={}): '.format(len(self.NSp['t']))+key_str
                description += '\n'+data_str
            if self.ENp:
                key_str = ', '.join(self.ENp.keys())
                data_str = 'ENpost (Nt={}): '.format(len(self.ENp['t']))+key_str
                description += '\n'+data_str
            if self.fieldstat:
                key_str = ', '.join(self.fieldstat.keys())
                data_str = 'BLfieldstat: '+key_str
                description += '\n'+data_str
        description += '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
        return description

class Grid(object):
    '''Three-dimensional grid used in SP-Wind'''
    
    def __init__(self,Nx2,Ny,Lx,Ly,zmeshf,Lref=1.0):

        self.__Nx2 = Nx2
        self.__Ny = Ny
        self.__Lx = Lx
        self.__Ly = Ly
        self.__xs = np.linspace(0.0,self.Lx,self.Nx2,endpoint=False)
        self.__ys = np.linspace(0.0,self.Ly,self.Ny,endpoint=False)

        Nz, zcc, zst = readsp.zmesh(zmeshf,Lref=Lref)
        self.__Nz = Nz
        self.__Lz = (zst[-1]-zst[0])
        self.__zcc = zcc
        self.__zst = zst

    @property
    def Nx2(self):
        return self.__Nx2

    @property
    def Ny(self):
        return self.__Ny

    @property
    def Nz(self):
        return self.__Nz

    @property
    def Lx(self):
        return self.__Lx

    @property
    def Ly(self):
        return self.__Ly

    @property
    def Lz(self):
        return self.__Lz

    @property
    def xs(self):
        return self.__xs

    @property
    def ys(self):
        return self.__ys

    @property
    def zcc(self):
        return self.__zcc

    @property
    def zst(self):
        return self.__zst

    def __str__(self):
        grid_str = '''\
            Numerical grid
            ****************************
            Domain size: {Lx:.1f} x {Ly:.1f} x {Lz:.1f}
            Grid size: {Nx2} x {Ny} x {Nz}'''.format(
            Lx=self.Lx,
            Ly=self.Ly,
            Lz=self.Lz,
            Nx2=self.Nx2,
            Ny=self.Ny,
            Nz=self.Nz)
        return textwrap.dedent(grid_str)

class ABL(object):
    '''ABL object containing all parameters'''

    def __init__(self,path=None,Lref=1.0,Uref=1.0,Tref=288.15):
        if not path:
            path = os.getcwd()

        tref = Lref/Uref

        nsfile = os.path.join(path,'NS.setup')
        enfile = os.path.join(path,'therm.setup')
        ekfile = os.path.join(path,'ekman.setup')
        ns_stp = readsp.ns_setup(nsfile,Lref=Lref,Uref=Uref)
        en_stp = readsp.en_setup(enfile,Lref=Lref,Uref=Uref)
        ek_stp = readsp.ek_setup(ekfile,Lref=Lref,Uref=Uref)
        self.__kappa = 0.4              #von Karman constant
        self.__z0 = ns_stp['z0']        #Surface roughness
        self.__dpdx = None              #Pressure gradient
        if ns_stp['force_mode']==1:
            self.__mode = 'MDBL'
        elif ns_stp['force_mode']==2:
            self.__mode = 'PDBL'
            self.__dpdx = ns_stp['force_param']
        
        self.__fc = None                #Coriolis parameter
        self.__gravity = None           #Gravitational acceleration
        self.__G = None                 #Geostrophic wind speed
        if ek_stp['ekman']:
            self.__fc = ek_stp['fc']
            self.__G  = ns_stp['force_param']
            self.__mode = 'EKM'
        if en_stp['thermo']:
            self.__gravity = en_stp['gravity']
            self.__mode = 'ABL'

    @property
    def mode(self):
        return self.__mode

    @property
    def kappa(self):
        return self.__kappa

    @property
    def z0(self):
        return self.__z0

    @property
    def dpdx(self):
        return self.__dpdx

    @property
    def fc(self):
        return self.__fc

    @property
    def gravity(self):
        return self.__gravity

    @property
    def G(self):
        return self.__G

    def __str__(self):
        abl_str='''\
            ABL parameters
            ****************************
            mode: {mode}
            kappa: {kappa}
            z0 = {z0:.4e} m'''.format(mode=self.mode,kappa=self.kappa,z0=self.z0)
        abl_str=textwrap.dedent(abl_str)
        if self.dpdx:
            abl_str=abl_str+'\n'+'dpdx = {dpdx:.3e} m/s^2'.format(dpdx=self.dpdx)
        if self.fc:
            abl_str=abl_str+'\n'+'fc = {fc:.4e} 1/s'.format(fc=self.fc)
        if self.gravity:
            abl_str=abl_str+'\n'+'g  = {gravity:.2f} m/s'.format(gravity=self.gravity)
        if self.G:
            abl_str=abl_str+'\n'+'G  = {G:.2f} m/s'.format(G=self.G) 
        return abl_str

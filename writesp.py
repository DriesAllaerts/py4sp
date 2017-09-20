#!/usr/bin/env python

'''
Module for writing setup files

Author: Dries Allaerts
Date: May 17, 2017
'''

import numpy as np
import os

def windfarm_setup(filename,Ct,layout,hubheight,radius,fwidth,Cp=0.0,tau=0.0,ALM_time_scale=100.0,ALM_hub_radius=0.015,ALM_rotational_speed=20.0,ALM_blade_pitch=0.0,Lref=1.0,Uref=1.0):
    Nt = layout.shape[0]
    layout = layout/Lref
    hubheight = hubheight/Lref
    radius = radius/Lref
    tau = tau*Uref/Lref
    with open(filename,'w') as file:
        file.write(str(Nt)+'\n')
        file.write('{:10.7f}'.format(Ct)+'{:15.7f}'.format(Cp)+'{:15.7f}'.format(tau)+'{:15.4f}'.format(ALM_time_scale)+'\n')
        for turb in range(Nt):
            file.write('{:10.7f}'.format(layout[turb,0])+'{:15.7f}'.format(layout[turb,1])+'{:15.7f}'.format(hubheight)+'{:15.7f}'.format(ALM_hub_radius)+'{:15.7f}'.format(radius)+'{:15.7f}'.format(ALM_rotational_speed)+'{:15.7f}'.format(ALM_blade_pitch)+'{:15.7f}'.format(fwidth)+'\n')

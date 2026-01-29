#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Methods to compute Heat Index (see Rothfusz et al. 1990).

"""

import numpy as np

def Heat_Index_simple(Tc, U):
    
    """
    T is air temperature (°C)
    U is relative humidity
    The coefficients are adapted to compute HI en °C
    
    """
    c1 = -8.78469476
    c2 = 1.61139411
    c3 = 2.33854884
    c4 = -0.14611605
    c5 = -1.230809e-2
    c6 = -1.642483e-2
    c7 = 2.21173e-3
    c8 = 7.2546e-4
    c9 = -3.58e-6
    HI = (c1+c2*Tc+c3*U+c4*Tc*U+c5*Tc**2+c6*U**2+c7*U*Tc**2+c8*Tc*U**2+c9*(U**2)*(Tc**2))
    
    return HI

def Heat_Index(Tc, U):
    """
    T is air temperature (°C)
    U is relative humidity
    The coefficients are adapted to compute HI en °C
    
    """
   
    HI=1.1*Tc+0.026*U-3.944
    
    HI=(HI+Tc)/2
    
    # new calculation where HI>=26.667°C
    
    ind_HI_80 = np.where(HI>=26.667)
            
    HI[ind_HI_80]=Heat_Index_simple(Tc[ind_HI_80],U[ind_HI_80])
    
    #adjustments
    
    ind_adjust_1=np.where((Tc<44.444) & (Tc>26.667) & (U<13))
    HI[ind_adjust_1]=HI[ind_adjust_1]-((13-U[ind_adjust_1])/7.2)*np.sqrt((17-np.abs(1.8*Tc[ind_adjust_1]-63.))/17)
    
    ind_adjust_2=np.where((Tc<30.556) & (Tc>26.667) & (U>85))
    HI[ind_adjust_2]=HI[ind_adjust_2]+((U[ind_adjust_2]-85)/18) * ((55-1.8*Tc[ind_adjust_2])/5)

    return HI

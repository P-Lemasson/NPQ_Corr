#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 10:42:14 2023

@author: P.Lemasson
"""


import numpy as np
import pandas as pd

def NPQ_Corr (depth, PAR, MLD, Chla , bbp):
    """
    Function to correct the effect of NPQ on the fluorescence signal 
    
    The method is based on the recommendations of Xing et al. 2018: adding a light threshold (PAR = 15) to existing methods such as Sackmann 2008.
    With the modification of Terrats 2020 consisting in applying the sigmoid proposed by Xing et al. 2018 for shallow mixing cases only under MLD.
    
    In deep mixing cases or in case of error with the sigmoid, the original method of Sackmann is applied, but above the light threshold.
    
    Chla must be smooth by rolling median window 11 (Xing 2018)
    PAR and bbp too
    Chla,bbp,PAR ==> same vertical resolution
    
    (Here sigmoide and Fratio are not smoothed)
    -----
    return 

    Fluorescence corrected for quenching
    Initial quenched fluorescence signal
    depth array
    ----
    
    """
    # Find Depth NPQ (ie zPAR = 15)
    zPAR15 = depth[np.absolute(PAR-15).argmin()]
    # MLD *0.9 : The factor of 0.9 ensures that possible biomass features near the pycnocline are avoided (Schallenberg 2022)
    MLD = MLD*0.9
    
    # Compute F-ratio from Sackmann 2008 and Sigmoide from Xing et al. 2018  (Only above zPAR15)
    Fratio = Chla[:np.absolute(PAR-15).argmin()+1]/bbp[:np.absolute(PAR-15).argmin()+1]
    
    
    XB18 = Chla[:np.absolute(PAR-15).argmin()+1] / (0.092 + 0.908 / (1 + (PAR[:np.absolute(PAR-15).argmin()+1] / 261)**2.2 ) )
    
    
    # =============================================================================
    # Shallow mixing cases ==> Under MLD sigmoide from Xing2018 and above correction from Sackmann2008
    # =============================================================================
    if zPAR15 > MLD :
        
        # Take Fratio-max as sigmoide/bbp at MLD 
        Fratio_max = XB18[np.absolute(depth - MLD).argmin()] / bbp[np.absolute(depth - MLD).argmin()]
    
        Fqc = np.concatenate( ( Fratio_max* bbp[: np.absolute(depth - MLD).argmin() +1] , 
                               XB18[np.absolute(depth - MLD).argmin() +1  :np.absolute(depth - zPAR15).argmin() +1 ]    ), axis=None )
        
        # =============================================================================
        # test if simoide method doesn't lead to smaller Fqc than Fq (that's impossible)
        # If it's the case for more than 80% : Sackmann method is applied
        # =============================================================================
        if len(np.where(Fqc < Chla[:len(Fqc)])[0]) / len(Fqc) > 0.8 :
            
            Fqc, depth_Fqc = SO8(depth , MLD, Fratio , bbp)
            
        else:
            
            depth_Fqc = depth[:len(Fqc)]
            
    else :
        
        Fqc, depth_Fqc = SO8(depth , MLD, Fratio , bbp)
    
    
    Fq = Chla
    
    return Fqc , depth_Fqc, Fq , Fratio, zPAR15


# seulement pour correction au dessus de MLD*0.9 ==> cas de deep mixing ou lorsque condition Catherine 
def SO8 (depth , MLD, Fratio , bbp) :
                
        
        Fratio_max = np.max(Fratio[:np.absolute(depth - MLD).argmin()+1])    # prend la plus haute valeur dans la couche d'eau au dessus de MLD
        
        depth_Fqc = depth[ : np.min(np.where( Fratio[:np.absolute(depth - MLD).argmin()+1] == Fratio_max))  +1 ]
        
        Fqc = Fratio_max * bbp[:len(depth_Fqc)]
               
        return  Fqc, depth_Fqc   









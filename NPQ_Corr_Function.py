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
    Function to Correct the Effect of NPQ on the Fluorescence Signal 
    
    This method is based on the recommendations of Xing et al. (2018), which involve adding a light threshold (PAR = 15) 
    to existing methods such as Sackmann (2008). 
    Additionally, it incorporates the modification proposed by Terrats (2020), which applies the sigmoid function suggested by Xing et al. (2018) 
    for shallow mixing cases only under the mixed layer depth (MLD).
    
    
    In cases of deep mixing or when encountering errors with the sigmoid function, the original method of Sackmann is applied.
    
    Chla, PAR, and bbp must undergo smoothing using a rolling median window of 11 (Xing, 2018).
    Chla,bbp,PAR ==> same vertical resolution
.
    (Note sigmoide and Fratio are not smoothed)
    -----
    return 

    Fluorescence corrected for quenching (Fqc)
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
        # test if sigmoide method doesn't lead to smaller Fqc than Fq (that's impossible)
        # If it's the case for more than 80% : Sackmann method is applied
        # =============================================================================
        if len(np.where(Fqc < Chla[:len(Fqc)])[0]) / len(Fqc) > 0.8 :
            # Application above MLD
            Fqc, depth_Fqc = SO8(depth , MLD, Fratio , bbp) 
            
        else:
            # Save depth array
            depth_Fqc = depth[:len(Fqc)]
            
    else :
        # Deep mixing case : SO8+ above zPAR15
        Fqc, depth_Fqc = SO8(depth , zPAR15 , Fratio , bbp) # Cas de deep mixing application de SO8+ donc au dessus de zPAR15
    
    
    Fq = Chla
    
    return Fqc , depth_Fqc, Fq , Fratio, zPAR15 , XB18



def SO8 (depth , profondeur_corr , Fratio , bbp) :
                
        # Find max value above MLD or zPAR=15
        Fratio_max = np.max(Fratio[:np.absolute(depth - profondeur_corr).argmin()+1])    
        # depth array
        depth_Fqc = depth[ : np.min(np.where( Fratio[:np.absolute(depth - profondeur_corr).argmin()+1] == Fratio_max))  +1 ]
        # Quenching correctiob
        Fqc = Fratio_max * bbp[:len(depth_Fqc)]
               
        return  Fqc, depth_Fqc   




def X12 (depth , profondeur_corr , Chla_q) : 
                
        # Extraction of the maximum CHLa concentration within the mixing layer (ML) 
        Chla_max = np.max(Chla_q [:np.absolute(depth - profondeur_corr).argmin()+1])
    
        depth_Fqc = depth[ : np.min(np.where( Chla_q[:np.absolute(depth - profondeur_corr).argmin()+1] == Chla_max) [0])  +1 ]
        
        # And extrapolation (or above zPAR15 in the case of X12+).
        Fqc = np.full(shape=len(depth_Fqc), fill_value=Chla_max)
               
        return  Fqc, depth_Fqc   




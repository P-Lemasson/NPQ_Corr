#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 12:16:55 2023

@author: oao2020


This script read Sprof files from floats 

Conditions : 
    
    Physique :
        Prise Adjusted si présent et si plus de 75% de QC = 1,2,5,8 (Si QC mauavais pas test du RT car sera forcément pas bon)
        Si Adjusted Absent : 
            Prise de RT et check si profil QC = A ou B
            
        Density is computes with the gsw module = TEOS10 & Read PAR data for conditons on Biogeo :
            If floats reach the MLD depth and the zPAR15
            or if PAR at the surface is higher than 15 (the NPQ limit from Xing 2018)
            (And if there is physics values at the surface)
    
    
    Biogeo : 
        Prise Adjusted si présent et si plus de 75% de QC = 1,2,5,8 (AU DELA DE zPAR15 = conseil de Louis car valeurs QC bon en dessous peuvent biaiser)
        Si Adjusted Absent : 
            Prise de RT si QC = 1 ou 2 Au dela de zPAR15 aussi
        CHL = Only raw data
        
        -Renosh Advise : Remove values QC 3 or 4 exspecially for PAR Adj (for compare with SOCA outputs) ------------ TO DO mais seulement sur le PAR
                        (No on physics variable because less important) 


"""


import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import lmfit
from lmfit.models import ExpressionModel

import toolbox_Argo as Aplt
import time_UTC_local as UTC

import gsw


#%% .nc read 
plt.close('all')

# =============================================================================
# Path to csv file wich contain path to .nc
# =============================================================================
path = '/home/oao2020/Documents/Pierrick/202212_BgcArgoSprof/dac/'
df = pd.read_csv(path + "SprofFiles_sorted.csv")

# List for save in a txt file at the end
alpha_NPQ_txt_X18_S08 = []
time_NPQ_txt_X18_S08 = []
latitude_txt_X18_S08 =[]
longitude_txt_X18_S08 = []
MLD_txt_X18_S08 = []
PAR_mean_txt_X18_S08 = []


# =============================================================================
# Loop : Read files BGC sorted .csv
# =============================================================================
# column = 'COR_CHL_DO_AUS'
column = 'CHL_DO'
#%% Conditions iterations
Data_Pres = 0
QC_Pres = 0

Data_SP = 0
QC_SP = 0

Data_Temp = 0
QC_Temp= 0

Data_PAR = 0
QC_PAR =0


Nuit = 0
Surface = 0
Visible = 0

Cond_zPAR = 0

Data_CHL= 0
Data_CHL_ADJ = 0
Len_CHL= 0
QC_CHL= 0

Data_Bbp = 0
Len_Bbp= 0
QC_Bbp= 0

Sigmoide = 0
fratio = 0

low_NPQ = 0
Len_NPQ = 0

Len_fit =0
soucis_fit = 0

R_cond = 0
NPQmax_positif = 0
slope_positive = 0
slope_maxvalue = 0

#%%
# random = np.random.randint(30,size=(1))
# for flotteur in range(len(df[column])):
for flotteur in [99]:
# for flotteur in random:
    data = xr.open_dataset(df[column][flotteur]) # prise des fichers du dac Coriolis avec fluo, radio et dans Autral
    
    # Extract some variables
    time=data['JULD'].values 
    lon = data['LONGITUDE'].values 
    lat = data['LATITUDE'].values 
    
    WMO = df[column][flotteur].split('/')[-1].split('_')[0]

    # List for save in a txt file at the end
    alpha_NPQ_X18_S08 = []
    MLD_NPQ_X18_S08 = []
    PAR_mean_X18_S08 = []
    time_NPQ_X18_S08 = []
    latitude_X18_S08 =[]
    longitude_X18_S08 = []
    
    # For load a random profile
    profile = np.random.randint(len(data['CYCLE_NUMBER'].values),size=(5))
  
# =============================================================================
# Loop by profile
# =============================================================================

    # for i in range(len(data['CYCLE_NUMBER'].values)): # loop for n profils
    # for i in profile:
    for i in [ 6]: #
        
#%%  PHYSIQUE Conditions
# =============================================================================
# Prise Adjusted si présent et si plus de 75% de QC = 1,2,5,8 (Si QC mauavais pas test du RT car sera forcément pas bon)
# Si Adjusted Absent : 
# Prise de RT et check si profil QC = A ou B      
# =============================================================================
        depth = data['PRES_ADJUSTED'].values[i]
        
        if not any(~np.isnan(depth)):
            
            depth = data['PRES'].values[i]
            
            if not any(~np.isnan(depth)):
                print(str(i) +'Pas de mesure de Pression Rt ni ADJ ' )
                Data_Pres += 1
                continue
            else:
            
                 depth_QC = data['PRES_QC'].values[i]
                 depth_QC = depth_QC[~pd.isnull(depth_QC)].astype('int')  
                 profile_QC = np.where ((depth_QC ==1) | (depth_QC==2)|(depth_QC ==5) | (depth_QC==8))
                 percent_depth = 100* len(profile_QC[0])/len(depth_QC)
                 
                 if not percent_depth > 75:
                   print(str(i) +': QC profil PRES trop mauvais : '+ str(percent_depth))
                   QC_Pres += 1
                   continue
        else : 
             depth_QC = data['PRES_ADJUSTED_QC'].values[i]
             depth_QC = depth_QC[~pd.isnull(depth_QC)].astype('int')  
             profile_QC = np.where ((depth_QC ==1) | (depth_QC==2)|(depth_QC ==5) | (depth_QC==8))
             percent_depth = 100* len(profile_QC[0])/len(depth_QC)
             
             if not percent_depth > 75:
               print(str(i) +': QC profil PRES_ADJ trop mauvais : '+ str(percent_depth))
               QC_Pres += 1
               continue
  
# =============================================================================
 
        SP = data['PSAL_ADJUSTED'].values[i]
        
        if not any(~np.isnan(SP)):
            
            SP = data['PSAL'].values[i]
            
            if not any(~np.isnan(SP)):
                print(str(i) +'Pas de mesure de PSAL Rt ni ADJ ' )
                Data_SP += 1
                continue
            
            else:
                SP_QC = data['PSAL_QC'].values[i]
                SP_QC = SP_QC[~pd.isnull(SP_QC)].astype('int')  
                profile_QC = np.where ((SP_QC ==1) | (SP_QC==2)|(SP_QC ==5) | (SP_QC==8))
                percent_SP = 100* len(profile_QC[0])/len(SP_QC)
                
                if not percent_SP > 75:
                  print(str(i) +': QC profil SP trop mauvais : '+ str(percent_SP))
                  QC_SP += 1
                  continue
                
        else : 
             SP_QC = data['PSAL_ADJUSTED_QC'].values[i]
             SP_QC = SP_QC[~pd.isnull(SP_QC)].astype('int')  
             profile_QC = np.where ((SP_QC ==1) | (SP_QC==2)|(SP_QC ==5) | (SP_QC==8))
             percent_SP = 100* len(profile_QC[0])/len(SP_QC)
             
             if not percent_SP > 75:
               print(str(i) +': QC profil SP ADJ trop mauvais : '+ str(percent_SP))
               QC_SP += 1
               continue
        
              
# =============================================================================      
        
        temp = data['TEMP_ADJUSTED'].values[i]
        
        if not any(~np.isnan(temp)):
            
            temp = data['TEMP'].values[i]
            
            if not any(~np.isnan(temp)):
                print(str(i) +'Pas de mesure de TEMP Rt ni ADJ ' )
                Data_Temp += 1
                continue
            else:
                Temp_QC = data['TEMP_QC'].values[i]
                Temp_QC = Temp_QC[~pd.isnull(Temp_QC)].astype('int')  
                profile_QC = np.where ((Temp_QC ==1) | (Temp_QC==2)|(Temp_QC ==5) | (Temp_QC==8))
                percent_Temp = 100* len(profile_QC[0])/len(Temp_QC)
                
                if not percent_Temp > 75:
                  print(str(i) +': QC profil TEMP trop mauvais : '+ str(percent_Temp))
                  QC_Temp += 1
                  continue
                
        else : 
             Temp_QC = data['TEMP_ADJUSTED_QC'].values[i]
             Temp_QC = Temp_QC[~pd.isnull(Temp_QC)].astype('int')  
             profile_QC = np.where ((Temp_QC ==1) | (Temp_QC==2)|(Temp_QC ==5) | (Temp_QC==8))
             percent_Temp = 100* len(profile_QC[0])/len(Temp_QC)
             
             if not percent_Temp > 75:
               print(str(i) +': QC profil TEMP trop mauvais : '+ str(percent_Temp))
               QC_Temp += 1
               continue
        

# =============================================================================
#    Compute density
# =============================================================================
    
        SA = gsw.SA_from_SP(SP,depth,lon[i],lat[i])
        CT = gsw.CT_from_t(SA,temp,depth)
        
        density = gsw.sigma0(SA,CT)
        density,depth_dens = Aplt.rem_nan(density,depth) # Soit profondeurs où il y a une mesure




#%% MLD / zPAR15

# =============================================================================
#  Read Available PAR for conditons
# =============================================================================
        PAR = data['DOWNWELLING_PAR_ADJUSTED'].values[i]# prise du PAR Adjusted
        PAR,depth2 = Aplt.rem_nan(PAR, depth[:len(PAR)])# liste PAR plus courte s'arrete plus tot
        PAR_data = '_ADJUSTED_QC'
        
        if not any(~np.isnan(PAR)):
            
            PAR = data['DOWNWELLING_PAR'].values[i]
            PAR,depth2 = Aplt.rem_nan(PAR, depth[:len(PAR)])
            PAR_data = '_QC'
            
            if not any(~np.isnan(PAR)):
                print(str(i) +': Pas de profil de PAR ni RT ni ADJ (descendant ?)')
                Data_PAR += 1
                continue
        

    
# =============================================================================
#  Nuit, densité en surface, profil atteint zPAR15 ou MLD        
# =============================================================================

        if not PAR[0]>15 :
            print(str(i) +': PAR de surface trop faible ie (<100) :'+str(PAR[0]) ) 
            Nuit += 1
            continue
        
        elif not any(depth_dens<10) or all(depth_dens<10):
            print(str(i) +' Pas de physique en surface : depth_dens[0] = ' + str(depth_dens[:0]))
            Surface += 1
            continue
        
 
        elif not density[-1]>density[np.max(np.where(depth_dens < 10)) +1 ]+0.03 or not any(PAR<15): 
            # Modif ici de or a and car probleme seulement si ne voit aucun des deux (permettrait plus de profils)
            # UPTADE : remise du or car si correction depuis MLD mais que PAR ne va pas profond on pourra pas avoir NPQsv
            # Mettre seulement condition de voir zPAR15 ? Non car besoin d'avoir la MLD en vue pour le PAR moyen
            print(str(i) +': Le flotteur ne voit pas la MLD ni le zPAR15')
            Visible += 1
            continue

        
# =============================================================================
# PAR process and zPAR15 determination                
# =============================================================================
   
    
    # Linear interpolation every 0.5m   
        PAR_inter, depth_inter_PAR  =Aplt.interpol(PAR, depth2, 0.5)
        
    # Find zPAR15   
        try :
            index_Zeu = np.min(np.where(PAR_inter < 15))
        except :
            try :
                index_Zeu = np.min(np.where(PAR_inter.astype(int) == 15))
            except:
                print('Prblm sur index zPAR15')
                Cond_zPAR += 1
                continue
            
        zPAR15 = depth_inter_PAR[index_Zeu]
        if zPAR15 < 5 :
            print('zPAR15 trop en surface : ' + str(zPAR15))
            Cond_zPAR += 1
            continue
            
    
        
#%% Chlorophylle

# =============================================================================
# Profils QC au dela de Zpar 15 (Pour les QC choix de prise celui de flotteur)
# ============================================================================= 
           
        CHLA_raw = data['CHLA'].values[i]
        CHLA_raw,depth1 = Aplt.rem_nan(CHLA_raw, depth)
         
        if not any(~np.isnan(CHLA_raw))  :
            print(str(i) +' Pas de profil de Fluo')
            Data_CHL += 1 
            continue
        elif len(CHLA_raw[:np.absolute(depth1 - zPAR15).argmin()])<3 :
             print(str(i) +'Fluo que <3 valeurs au dessus de zPAR15')
             Len_CHL += 1 
             continue
         
        CHL_QC =   data['CHLA_QC' ].values[i] 
        CHL_QC = CHL_QC[~pd.isnull(CHL_QC)].astype('int')  
        # Pour la CHL on prend aussi le flag 3 car représente valeurs ajustables" et donc celles quenchées
        profile_QC = np.where((CHL_QC[:np.absolute(depth1 - zPAR15).argmin()]   ==0) |(CHL_QC[:np.absolute(depth1 - zPAR15).argmin()]   ==1) | 
                                (CHL_QC[:np.absolute(depth1 - zPAR15).argmin()]==2) | (CHL_QC[:np.absolute(depth1 - zPAR15).argmin()]   ==3)) 
        percent_CHL = 100* len(profile_QC[0])/len(CHL_QC[:np.absolute(depth1 - zPAR15).argmin()])
        
        if not percent_CHL > 75:
           print(str(i) +': QC profil CHL trop mauvais : '+ str(percent_CHL))
           QC_CHL+= 1
           continue
    
        
  
  
 # =============================================================================
 # Application des corrections dark et slope sur Chlraw 
 # =============================================================================
        # CHLA_adj = data['CHLA_ADJUSTED'].values[i] # Read Adjusted CHL
        # CHLA_adj = CHLA_adj[~np.isnan(CHLA_adj)]  # Remove all NaN values
 
        # try :
        #     index = np.min(np.where(CHLA_adj != CHLA_adj[0])) # Take the index when below the NPQ correction for find the dark corr
        # except:
        #     print('Prblm avec CHLA Adjusted')
        #     Data_CHL_ADJ += 1
        #     continue
       
        # CHLA_raw = CHLA_raw*0.5# Slope correction
        # # Dark correction
        # dark_corr =  CHLA_raw[index]-CHLA_adj[index]
        # if dark_corr>0 :
        #     CHLA_raw = CHLA_raw - abs(dark_corr)
        # elif dark_corr<0 :
        #       CHLA_raw = CHLA_raw + abs(dark_corr)
         
                 
 # =============================================================================
 #  Interpolation and smooth
 # =============================================================================
 
 
         # Linear interpolation every 0.5m    
        Chlraw_inter, depth_inter_CHL  =Aplt.interpol(CHLA_raw, depth1, 0.5)
        
        # Filtre des données se situant au dessus de zPAR15 (fait dans la fonction Ryan Keogh pour valeurs nulles et negatives)
        # => 0.02 limite de detection si plus faible (vraiment utile au finale je n'y met pas dans le rapport car cette limite a du sens sans les différentes correction de l'ADJ ?)
        Chlraw_inter[np.where((Chlraw_inter[:index_Zeu] <= 0.02) & (Chlraw_inter[:index_Zeu]>0)) ] = 0.02    
        
        
        #Smooth mobil median 11 => Advise from XING 2018
        Chlraw_smth_inter = Aplt.median_smooth(Chlraw_inter, 11) 
     
        # Nan values at surface => I assume that it's constant in the first 2 meter 
        index = np.where(np.isnan(Chlraw_smth_inter[:10]))
        Chlraw_smth_inter[index] = Chlraw_smth_inter[np.max(index)+1] 
     
        # Nan values at surface => Interpolate values
        # Chlraw_smth_inter[np.isnan(Chlraw_smth_inter)] = Chlraw_inter[np.isnan(Chlraw_smth_inter)]
 

 # =============================================================================
 # Vision du CHLA_raw smth      
 # =============================================================================
     # plt.figure()
     # plt.plot(CHLA_raw,depth1)
     # # plt.plot(Chlraw_smth_inter,depth_inter_CHL_smth_CHL)
     # # plt.plot(Chlraw_inter,depth_inter_CHL_CHL)
     # plt.gca().invert_yaxis()
     
        
    
 
    
# ============================================================================= 
#%% Bbp  


# =============================================================================
# Profils QC au dela de Zpar 15 
# =============================================================================  
  
        Bbp = data['BBP700_ADJUSTED'].values[i]
        Bbp,depth3 = Aplt.rem_nan(Bbp, depth)
        Bbp_data = '_ADJUSTED_QC'
        
        if not any(~np.isnan(Bbp)) or len(Bbp[:np.absolute(depth3 - zPAR15).argmin()])<3 : #Si moins de 3 valeurs au dessus de zPAR15

            Bbp = data['BBP700'].values[i]
            Bbp,depth3 = Aplt.rem_nan(Bbp, depth)
            Bbp_data = '_QC'
            
            if not any(~np.isnan(Bbp))  :
                Data_Bbp += 1 
                print(str(i) +' Pas de profil de Bbp ')
                continue
            elif len(Bbp[:np.absolute(depth3 - zPAR15).argmin()])<3 :
                print(str(i) + 'que <3 valeurs (ADJ et RT)')
                Len_Bbp+= 1
                continue

        Bbp_QC =   data['BBP700'+ Bbp_data ].values[i] #Prise des QC du profils adj ou Rt selon choix precedent
        Bbp_QC = Bbp_QC[~pd.isnull(Bbp_QC)].astype('int')  
        
        if Bbp_data == '_ADJUSTED_QC':
            profile_QC = np.where((Bbp_QC[:np.absolute(depth3 - zPAR15).argmin()]   ==1) | (Bbp_QC[:np.absolute(depth3 - zPAR15).argmin()]==2)| 
                                  (Bbp_QC[:np.absolute(depth3 - zPAR15).argmin()]==5)| (Bbp_QC[:np.absolute(depth3 - zPAR15).argmin()]==8)) #prise de surface à zPAR15 | = or
            percent_Bbp = 100* len(profile_QC[0])/len(Bbp_QC[:np.absolute(depth3 - zPAR15).argmin()])
            
        else:
            profile_QC = np.where((Bbp_QC[:np.absolute(depth3 - zPAR15).argmin()]   ==0) |(Bbp_QC[:np.absolute(depth3 - zPAR15).argmin()]   ==1) | 
                                  (Bbp_QC[:np.absolute(depth3 - zPAR15).argmin()]==2)) #prise de surface à zPAR15 | = or
            percent_Bbp = 100* len(profile_QC[0])/len(Bbp_QC[:np.absolute(depth3 - zPAR15).argmin()])
        
        if not percent_Bbp > 75:
           print(str(i) +': QC profil Bbp trop mauvais : '+ str(percent_Bbp))
           QC_Bbp+= 1
           continue

# =============================================================================
#  Interpolation and smooth
# =============================================================================
   
        # Linear interpolation every 0.5m
        Bbp_inter, depth_inter_Bbp  =Aplt.interpol(Bbp, depth3, 0.5) 
        
        #Smooth mobil median 11 
        Bbp_smth_inter = Aplt.median_smooth(Bbp_inter, 11)
                
        # Nan values at surface => I assume that it's constant in the first 2 meter 
        index = np.where(np.isnan(Bbp_smth_inter[:10]))
        Bbp_smth_inter[index] = Bbp_smth_inter[np.max(index)+1] 
        
        # Nan values at surface => Interpolate values
        # Bbp_smth_inter[np.isnan(Bbp_smth_inter)] = Bbp_inter[np.isnan(Bbp_smth_inter)]
            
            
# =============================================================================
# Vision du Bbp       
# =============================================================================
        # plt.figure()
        # plt.plot(Bbp,depth3)
        # plt.plot(Bbp_inter,depth_inter_Bbp)
        # plt.plot(Bbp_smth_inter,depth_inter_Bbp)
            
        # # Bbp_adj = data['BBP700_ADJUSTED'].values[i]
        # # Bbp_adj,depth4 = Aplt.rem_nan(Bbp_adj, depth)
        # # plt.plot(Bbp_adj,depth4)
       
        # x1,x2,y1,y2 = plt.axis()  
        # plt.axis((x1,x2,0,200))
        
        # plt.gca().invert_yaxis()



# =============================================================================
#%% PAR
        PAR_QC =   data[ 'DOWNWELLING_PAR' + PAR_data ].values[i] #Prise des QC du profils adj ou Rt selon choix precedent
        PAR_QC = PAR_QC[~pd.isnull(PAR_QC)].astype('int')  
        
        if PAR_data == '_ADJUSTED_QC':
            profile_QC = np.where((PAR_QC[:np.absolute(depth2 - zPAR15).argmin()]   ==1) | (PAR_QC[:np.absolute(depth2 - zPAR15).argmin()]==2)|
                                  (PAR_QC[:np.absolute(depth2 - zPAR15).argmin()]==5)| (PAR_QC[:np.absolute(depth2 - zPAR15).argmin()]==8)) #prise de surface à zPAR15 | = or
            percent_PAR = 100* len(profile_QC[0])/len(PAR_QC[:np.absolute(depth2 - zPAR15).argmin()])
            
            # Remove PAR values with QC = 3 or 4 
            PAR = np.delete(PAR , np.where((PAR_QC==3) & (PAR_QC == 4)), axis=None    )
            depth2 = np.delete(depth2 , np.where((PAR_QC==3) & (PAR_QC == 4)), axis=None    )
            
            
        else : 
            profile_QC = np.where((PAR_QC[:np.absolute(depth2 - zPAR15).argmin()]   ==0) |(PAR_QC[:np.absolute(depth2 - zPAR15).argmin()]   ==1) | 
                                  (PAR_QC[:np.absolute(depth2 - zPAR15).argmin()]==2))
            percent_PAR = 100* len(profile_QC[0])/len(PAR_QC[:np.absolute(depth2 - zPAR15).argmin()])
        
            # Remove PAR values with QC = 4 
            PAR = np.delete(PAR , np.where(PAR_QC == 4), axis=None    )
            depth2 = np.delete(depth2 , np.where(PAR_QC == 4), axis=None    )
        
        if not percent_PAR > 75:
           print(str(i) +': QC profil PAR trop mauvais : '+ str(percent_PAR))
           QC_PAR += 1
           continue
        
        else :
            print(i)

            
             # Linear interpolation every 0.5m of the new PAR without QC=4 (or 3)
            PAR_inter, depth_inter_PAR  =Aplt.interpol(PAR, depth2, 0.5)
            #Smooth mobil median 11 
            PAR_smth_inter = Aplt.median_smooth(PAR_inter, 11) 
            
            #Remplace nan value by the inter (Keep 'natural' value at surface)
            PAR_smth_inter[np.isnan(PAR_smth_inter)] = PAR_inter[np.isnan(PAR_smth_inter)] 
            
            
          
             
            
            
        
# =============================================================================
# Vision du PAR       
# =============================================================================
            # plt.figure()
            # # plt.plot(PAR,depth2,'--')
            # plt.plot(PAR_smth_inter,depth_inter_PAR)
            # plt.plot(PAR_inter,depth_inter_PAR)
            
            # # PAR_adj = data['PAR_ADJUSTED'].values[i]
            # # PAR_adj,depth4 = Aplt.rem_nan(PAR_adj, depth)
            
            # plt.axhline(y= zPAR15,color='red',label='Z PAR=15')
            # plt.gca().invert_yaxis()
# =============================================================================        
       
      
       
        



#%% MLD TEOS10 AND mean PAR     
 
            MLD = Aplt.MLD(density, depth_dens)
           
# =============================================================================
# Vision de la MLD      
# =============================================================================
            # plt.figure()
            # plt.plot(density,depth_dens)
            # plt.axhline(y= MLD,color='grey',label='MLD')
            # plt.gca().invert_yaxis()
            
# =============================================================================
#  Le PAR moyen dans la MLD  
# =============================================================================
            index = np.absolute(depth_inter_PAR - MLD).argmin()
            
            PAR_inter[np.where(PAR_inter<0)[0]] = 0
            
            if MLD > max(depth_inter_PAR) and abs(MLD - max(depth_inter_PAR)) > 5 : # Si MLD plus profonde que prof max PAR de 5metre au moins (sinon ballec)
            
                if min(PAR_inter) < PAR_inter[0]/100 : #Si plus petite valeur du PAR est plus faible que 1% de la surface
                
                # Attribution de 0 en desous jusqu'a MLD
                    Nbr_add_values = int((MLD - max(depth_inter_PAR))/0.5)
                     
                    PAR_inter_median_calcul = np.hstack ( [ PAR_inter  , np.linspace(0, 0 , Nbr_add_values)]   ) # permet de faire extend mais avec ndarray
                
                    PAR_moyen = np.median(PAR_inter_median_calcul) # mediane de tout du coup vu que deeper value = a MLD
                    
                else:
                    PAR_moyen = np.nan # Attribution de Nan car ne va pas jusqu'a MLD et n'atteint pas 1% de PAR surface et donc pas moyen attribuer des 0
            
            
            else :
                PAR_moyen = np.median(PAR_inter[:index+1])         # Pas de prblm les mesures de PAR voient MLD    

        

            
#%% Sigmoide 
            # Sigmoide from Xing2018 for Shallow mixing cases
            try:
                XB18 = Chlraw_inter[:len(PAR_smth_inter[:index_Zeu+1])] / (0.092 + 0.908 / (1 + (PAR_smth_inter[:index_Zeu+1] / 261)**2.2 ) ) # On prend index Zeu car on veut eviter les valeurs negatives de PAR en dessous, numpy n'aime pas trop les valeur negative puissance de float
            except:
                print('Probleme Sigmoide')
                Sigmoide+= 1
                continue
            
            # plt.figure()
            # plt.plot(XB18 ,depth_inter_PAR )
            
            
#%% Fratio
            # F ratio (Sackmann 2008) = CHL / Bbp 
            try :
                Fratio = Chlraw_smth_inter[:len(PAR_smth_inter[:index_Zeu+1])] / Bbp_smth_inter[:len(PAR_smth_inter[:index_Zeu+1])] 
            except:
                print('Fratio soucis')
                fratio += 1
                continue
            # For NPQ correction smooth advise from Louis
            Fratio_smth = Aplt.median_smooth(Fratio, 11) 
            
            #Remplace nan value by the inter (Keep 'natural' value at surface)
            Fratio_smth[np.isnan(Fratio_smth)] = Fratio[np.isnan(Fratio_smth)]
            


#%% NPQ correction Sackmann 2008

# =============================================================================
#   Application de méthode de Sackmann 2008 (avec le Bbp) + zPAR15 = SO8+
#   This method originating from Sackmann et al., (2008) it assumes uniformity in particle composition inside the mixed layer
# =============================================================================
            def SO8 (zPAR15 , MLD, depth_inter_Bbp , Fratio_smth) :
                
                        if  zPAR15 > MLD*0.9 : # correction depuis MLD*0.9 (shallow mixing case)
                            zCORR = MLD*0.9
                        else:
                            zCORR = zPAR15 
                            
                        # calculate the difference array # find the index of minimum element from the array
                        index = np.absolute(depth_inter_Bbp - zCORR).argmin()
    
                        Fratio_max = np.max(Fratio_smth[:index+1])    # prend la plus haute valeur dans la couche d'eau
                        
                        index = np.min(np.where(Fratio_smth[:index+1] == Fratio_max))
                        depth_NPQ_corr = depth_inter_Bbp[index]
                        
                               
                        return  Fratio_max,  depth_NPQ_corr # Valeur max et sa profondeur
                    
                    
            Fratio_max, depth_NPQ_corr = SO8(zPAR15, MLD, depth_inter_Bbp[:len(Chlraw_smth_inter)], Fratio_smth) 
            
            Fratio_depth_NPQ_corr = depth_inter_CHL[: np.min(np.where(depth_inter_CHL== depth_NPQ_corr))+1 ] # reconstruction depth array
            Fratio_NPQ_corr = Fratio_max * Bbp_smth_inter[:len(Fratio_depth_NPQ_corr)] # reconstruction Fratio array

#%% NPQ correction X18_SO8+ from Terrats 2020

            def X18_S08 (zPAR15 , MLD, depth_inter_Bbp  , XB18 , Bbp_smth_inter) :
                                    
                        if  zPAR15 > MLD*0.9 : #  (shallow mixing case) Under MLD*09 Application : XB18 and above Sackamnn 2008
                        
                            
                            index = np.absolute(depth_inter_Bbp - MLD*0.9).argmin() # index of MLD*0.9 position
                            
                            Fratio_max = XB18[index] / Bbp_smth_inter[index] # Take XB18 when cross MLD*0.9 for correction
                            
                            # Depth array
                            X18_S08_depth_NPQ_corr = depth_inter_Bbp[: np.absolute(depth_inter_Bbp - zPAR15).argmin() +1  ]          
                            # Corrected array : Fratio max * bbp above MLD 09 and XB18 below
                            X18_S08_NPQ_corr = np.concatenate( ( Fratio_max* Bbp_smth_inter[: np.absolute(depth_inter_Bbp - MLD*0.9).argmin() +1] , 
                                                XB18[np.absolute(depth_inter_Bbp - MLD*0.9).argmin() +1  :np.absolute(depth_inter_Bbp - zPAR15).argmin() +1 ]    ), axis=None )
                            
                            return  X18_S08_NPQ_corr,  X18_S08_depth_NPQ_corr
                        
                        
                        else:  #Correction above zPAR15 : Sackmann 2008
                    
                    
                            index = np.absolute(depth_inter_Bbp - zPAR15).argmin()
                           
                            Fratio_max = np.max(Fratio_smth[:index+1])    # Max value of F-ratio above zPAR15
                            
                            index = np.min(np.where(Fratio_smth[:index+1] == Fratio_max))  # find depth of max F ratio
                            depth_NPQ_corr = depth_inter_Bbp[index]
                            
                            # Depth array
                            X18_S08_depth_NPQ_corr = depth_inter_Bbp[: np.min(np.where(depth_inter_Bbp== depth_NPQ_corr))+1 ]
                            # Corrected array
                            X18_S08_NPQ_corr = Fratio_max * Bbp_smth_inter[:len(X18_S08_depth_NPQ_corr)]
                        
                            return X18_S08_NPQ_corr,  X18_S08_depth_NPQ_corr
                    

            X18_S08_NPQ_corr,  X18_S08_depth_NPQ_corr  = X18_S08(zPAR15, MLD, depth_inter_PAR, XB18, Bbp_smth_inter)



#%% Compute NPQsv
            def Ryan_Keogh2020 (CHL_corr , CHL , PAR)  :
                Fqc = np.copy(CHL_corr) 
                Fq = np.copy(CHL[:len(CHL_corr)]) 
                
                # on delete valeur CHL_inter si =0 ou <0
                index = np.where(Fq <= 0)
                Fqc = np.delete(Fqc,index)
                Fq = np.delete(Fq,index)
                
                NPQsv = (Fqc - Fq) /Fq
                PARsv = np.delete(np.copy(PAR[:len(CHL_corr)]),index)
               
                
                return NPQsv , PARsv
    
            # NPQsv_SO8 , PARsv_SO8 = Ryan_Keogh2020( Fratio_NPQ_corr , Chlraw_inter, PAR_inter)
            # NPQsv_X18_S08, PARsv_X18_S08 = Ryan_Keogh2020( X18_S08_NPQ_corr , Chlraw_inter, PAR_inter)
    
            # Use CHL an PAR smth for compare with SOCA outputs
            NPQsv_SO8 , PARsv_SO8 = Ryan_Keogh2020( Fratio_NPQ_corr , Chlraw_smth_inter, PAR_smth_inter)
            NPQsv_X18_S08, PARsv_X18_S08 = Ryan_Keogh2020( X18_S08_NPQ_corr , Chlraw_smth_inter, PAR_smth_inter)
            

            

#%% Conditions if from Schmetig 

# S'applique si la correction est inférieur au signal raw, ce qui n'est pas possible
# Arrive lorsque correction comme louis avec sigmoide en dessous de MLD (donc shallow mixing)           
# Que le Bbp est relativement Constant au dessus de zPAR15 et que la fluo raw augmente au dessus de zPAR15 
# = un cas très spécifique
          
# Donc si 80% de Fqc > Fq application de SO8+  
            if zPAR15 > MLD*0.9 and  len(np.where(X18_S08_NPQ_corr < Chlraw_inter[:len(X18_S08_NPQ_corr)])[0]) / len(X18_S08_NPQ_corr) > 0.8 :
                # Si existe on l'affiche
                # plt.figure(figsize=[12,8])
                
                
                # plt.plot(PAR_inter,depth_inter_PAR,'red',label='PAR_Res0.5',alpha= 0.25)
                # # plt.plot(Bbp_smth_inter,depth_inter_Bbp,'red',label='PAR_inter0.5')
                
                # plt.axis([0,np.max(PAR_inter),0,200])
                # plt.xlabel('PAR ',loc='center',fontweight='semibold')
                
                # plt.legend(loc='upper right')
                
                
                # plt.twiny()
                
                
                # plt.plot(CHLA_raw,depth1,'.', label='CHLA_raw') 
                
                # # plt.plot(Chlraw_smth_inter,depth_inter_CHL,color='darkorange', label='CHLA_smth17_Resolution0.5')
                # plt.plot(Chlraw_smth_inter,depth_inter_CHL,color='green', label='CHLA_smth_Res0.5')
                
                # # plt.plot(CHL_NPQ_corr_inter,depth_NPQ_corr_inter,color='darkorange',linestyle='dashed', label='X12+ NPQ correction')
                
                # # plt.plot(Fratio_NPQ_corr   ,Fratio_depth_NPQ_corr,linestyle='dashed',color='k',label = 'SO8+ NPQ correction')
                
                # # if  zPAR15 > MLD*0.9 :
                # plt.plot(XB18 ,depth_inter_PAR[:index_Zeu+1],'bo',markersize = 2.5, alpha= 0.5,label='Sigmoide X18')
                # plt.plot(X18_S08_NPQ_corr,X18_S08_depth_NPQ_corr,'cornflowerblue',label='X18_SO8 NPQ correction',linestyle='dashed')
                
                # plt.xlabel('CHL [mg/m3]',loc='center',fontweight='semibold')
                # plt.ylabel('Pressure (dbar)',fontweight='semibold')
                # plt.legend(loc='lower right')
                
                # plt.axhline(y= zPAR15,color='red',label='Z PAR=15')
                
                # plt.annotate('zPAR = 15',
                # xy=(0, zPAR15), xycoords='data',
                # xytext=(0, zPAR15-5), fontsize=8,color='red')
                
                
                # plt.axhline(y= MLD,color='black')
                # plt.axhline(y= MLD*0.9,color='grey',label='MLD*0.9')
                
                # plt.annotate('MLD*0.9',
                #       xy=(0, MLD*0.9), xycoords='data',
                #       xytext=(0, MLD*0.9-5), fontsize=8)
                
                # plt.annotate('MLD',
                #       xy=(0, MLD*0.9), xycoords='data',
                #       xytext=(0, MLD+8), fontsize=8) 
                
                # plt.axis([0,np.max([np.max(Chlraw_inter),np.max(Fratio_NPQ_corr),  np.max(X18_S08_NPQ_corr) ]  ),0,200])
                
                
                
                # plt.gca().invert_yaxis()
                # # plt.xlabel('PAR [microMoleQuanta/m^2/sec]',loc='center',fontweight='semibold')
                
    
                
                # # plt.suptitle ('Profil vertical fluo {0}'.format(data['CYCLE_NUMBER'].values.astype(int)[i]) +' {0} '.format(dire)+WMO + '  '+ d)
                
                
                
                
                
                
                NPQsv_X18_S08 = NPQsv_SO8
                PARsv_X18_S08 = PARsv_SO8
                print('Fqc avec X18_SO8 mauvais application de So8 ..................................................')




#%% Alpha fit
            if all(NPQsv_X18_S08 < 0.10) :
                print(str(i) +' Pas de NPQsv>0.1 (Correction non necessaire ?) ' ) 
                low_NPQ += 1
                continue
            elif len(NPQsv_X18_S08)<10 :
                print(str(i) +' len NPQ <10 (Correction non necessaire ?) ' ) 
                Len_NPQ += 1
                continue
# # =============================================================================
# #    Pour plot les points sans l'offset de Shallenberg 
# # =============================================================================
           
            NPQsv_plot , PARsv_plot = NPQsv_X18_S08 , PARsv_X18_S08
            # On leur vire quand meme les valeurs negatives = juste visuel
            index = np.where(NPQsv_plot<0)
            NPQsv_plot = np.delete(NPQsv_plot,index)
            PARsv_plot = np.delete(PARsv_plot,index)
            
            depth_NPQ = np.delete(X18_S08_depth_NPQ_corr,index) # For visualize NPQ with depth
# # =============================================================================
# #     Shallenberg : offset pour aller vers ordonnée origine
# # =============================================================================
            def Shall_offset (NPQsv,PARsv):
                
                # on supp lorsque NPQsv<0.1 (aussi ceux negatifs par meme occasion)
                index = np.where(NPQsv<0.1)
                NPQsv = np.delete(NPQsv,index)
                PARsv = np.delete(PARsv,index) 
                
                # offset la valeur au niveau de NPQ = 0.1 (avoir ordonnée origine après offset)
                NPQsv_offset = NPQsv[-1]
                PARsv_offset = PARsv[-1]
                NPQsv = NPQsv - NPQsv[-1]
                PARsv = PARsv -PARsv[-1] 
                
                #rm negative value , happen when [-1] > à d'autre valeurs plus au dessus
                PARsv= PARsv[NPQsv>=0]
                NPQsv = NPQsv[NPQsv>=0]
                
                NPQsv= NPQsv[PARsv>=0]
                PARsv = PARsv[PARsv>=0]
                
                # # # Pour le fit avoir les points croissant (depart effet visuel) mais affecte le NPQmax (seulement lui car grosse variation du PAR en surface et non pas plus en profondeur où alpha interet)
                idx = np.argsort(PARsv)
                NPQsv = NPQsv[idx]
                PARsv = PARsv[idx]
                
                return NPQsv,PARsv, NPQsv_offset, PARsv_offset
              

            NPQsv_X18_S08 , PARsv_X18_S08,  NPQsv_offset_X18_S08, PARsv_offset_X18_S08 = Shall_offset(NPQsv_X18_S08 , PARsv_X18_S08)
            
            
            
            
#%% Fit Platt 80             

            from scipy.optimize import curve_fit # = leats square method
            
            def Platt(x,NPQmax ,slope  ):
                y = NPQmax * (1-np.exp(-slope*x/NPQmax)) 
                return y
            
# =============================================================================   
            if len(NPQsv_X18_S08) < 10 :
                print('Pas assez de valeur pour un fit')
                Len_fit += 1
                continue

# =============================================================================         
            try :
                parameters, covariance = curve_fit(Platt, PARsv_X18_S08, NPQsv_X18_S08, p0 = [3,0.01])  #p0 = initial guess default = 1
            except:
                print(str(i) + 'Soucis avec le fit')
                soucis_fit += 1
                continue
            NPQmax = parameters[0]
            slope = parameters[1]
            fit = Platt(PARsv_X18_S08, NPQmax, slope)
            
            
            R2 = 1 -  (sum( (NPQsv_X18_S08- fit) **2)  /  sum( (NPQsv_X18_S08 - np.mean(NPQsv_X18_S08)) **2)  )  #le R squared
            
            
            
            MSE = sum( (NPQsv_X18_S08- fit) **2) / len(NPQsv_X18_S08)
            RMSE = MSE **0.5
# =============================================================================
           
            
            t = UTC.time_UTCtolocal(time[i],lon[i]) # UTC to local hour
            ts = pd.to_datetime(str(t)) 
            d = ts.strftime('%Y/%m/%d %H:%M:%S') # make the date easier to read
            dire = data['DIRECTION'].values[i]
            
            if R2>=0.8 and NPQmax >0 and slope <= 0.1 and slope>0:  #and slope<0.1 vraiment utile ? oui  il l'est (voir screen) = Condition comme Shallenberg
      
                alpha_NPQ_X18_S08.append(round(slope,4))
                time_NPQ_X18_S08.append( ts.strftime('%Y/%m/%d')) 
                
                latitude_X18_S08.append(lat[i])
                longitude_X18_S08.append(lon[i]) 

                MLD_NPQ_X18_S08.append(MLD)
                PAR_mean_X18_S08.append(PAR_moyen)
                
            elif R2 < 0.8:
                R_cond+= 1
            elif NPQmax <=0 :
                NPQmax_positif += 1
            elif slope > 0.1:
                 slope_maxvalue+= 1
            elif slope<= 0:
                 slope_positive += 1
#%% Figures 

            plt.figure(figsize=[12,8])
            plt.subplot(131)
            
            # plt.plot(PAR_inter,depth_inter_PAR,'red',label='PAR_Res0.5',alpha= 0.5)
            plt.plot(PAR_smth_inter,depth_inter_PAR,'red',label='PAR_smth_Res0.5',alpha= 0.5)
            
            plt.axis([0,np.max(PAR_inter),0,200])
            plt.xlabel('PAR ',loc='center',fontweight='semibold')
            plt.ylabel('Pressure (dbar)',fontweight='semibold')
            plt.legend(loc='lower left')
            
            
            plt.twiny()
            
            
            plt.plot(CHLA_raw,depth1,'go', label='CHLA_raw',markersize=5,alpha=0.5) 
            
            # plt.plot(Chlraw_smth_inter,depth_inter_smth_CHL,color='darkorange', label='CHLA_smth17_Resolution0.5')
            plt.plot(Chlraw_inter,depth_inter_CHL,color='green', label='CHLA_Res0.5')
            
            
            
            # plt.plot(XB18 ,depth_inter_PAR[:index_Zeu+20],'bo',markersize = 2.5, alpha= 0.5,label='Sigmoide X18')
            plt.plot(X18_S08_NPQ_corr,X18_S08_depth_NPQ_corr,'cornflowerblue',label='X18_SO8 NPQ correction',linestyle='dashed')
            
            plt.xlabel('CHL [mg/m3]',loc='center',fontweight='semibold')
            plt.ylabel('Pressure (dbar)',fontweight='semibold')
            plt.legend(loc='lower right')
            
            plt.axhline(y= zPAR15,color='red',label='Z PAR=15')
            plt.annotate('zPAR = 15',xy=(0, zPAR15), xycoords='data',xytext=(0.05, zPAR15-5), fontsize=8,color='red')
            
            
            plt.axhline(y= MLD,color='black')
            plt.axhline(y= MLD*0.9,color='grey',label='MLD*0.9')
            
            plt.annotate('MLD*0.9',xy=(0, MLD*0.9), xycoords='data', xytext=(0.05, MLD*0.9-5), fontsize=8)
            
            plt.annotate('MLD',xy=(0, MLD*0.9), xycoords='data', xytext=(0.05, MLD+8), fontsize=8) 
            
            plt.axis([0,np.max([np.max(Chlraw_inter),np.max(Fratio_NPQ_corr),  np.max(X18_S08_NPQ_corr) ]  ),0,150])
            
            
            
            plt.gca().invert_yaxis()
            
            
            
            
# =============================================================================  # =============================================================================   
            plt.subplot(132)


            plt.plot(Bbp_smth_inter,depth_inter_Bbp,'k',label='bbp_smth_res0.5',alpha= 0.5)
            
            
            plt.xlabel('Bbp ',loc='center',fontweight='semibold')
           
            plt.legend(loc='lower left')
            
            
            
            
            plt.twiny()
            
            plt.plot(Chlraw_smth_inter,depth_inter_CHL,color='darkorange', label='CHLA_smth11_Res0.5')
            plt.plot(X18_S08_NPQ_corr,X18_S08_depth_NPQ_corr,'cornflowerblue',label='X18_SO8 NPQ correction',linestyle='dashed')
            
            
            plt.axhline(y= zPAR15,color='red')
            plt.annotate('zPAR = 15',xy=(0, zPAR15), xycoords='data',xytext=(0.05, zPAR15-5), fontsize=8,color='red')
            
            plt.axhline(y= MLD*0.9,color='grey')
            plt.annotate('MLD*0.9',xy=(0, MLD*0.9), xycoords='data',xytext=(0.05, MLD*0.9-5), fontsize=8)
            
             
            plt.axis([0,np.max([np.max(Chlraw_inter),np.max(Fratio_NPQ_corr),  np.max(X18_S08_NPQ_corr) ]  ),0,150])
            plt.gca().invert_yaxis()


            plt.xlabel('CHL [mg/m3]',loc='center',fontweight='semibold')
            
            plt.legend(loc='lower right')
# =============================================================================  # =============================================================================


            plt.subplot(133)
            
          
            plt.plot(PARsv_plot,NPQsv_plot,'o',color='cornflowerblue',alpha=0.5)
            
            plt.axis([0,np.max(PARsv_X18_S08   )+PARsv_offset_X18_S08 ,0, np.max(NPQsv_X18_S08   )+NPQsv_offset_X18_S08])
            plt.xlabel('iPAR [microMoleQuanta/m^2/sec]',fontweight='semibold')
            plt.ylabel('NPQsv',fontweight='semibold')




            plt.plot(PARsv_X18_S08+PARsv_offset_X18_S08, fit +NPQsv_offset_X18_S08, 'k--', label='X18_S08_fit')
            
            plt.annotate('R2 =' + str(round(R2,4))+' \n' + 'RMSE =' + str(round(RMSE,4)) +' \n' + r'$\alpha$ =' +str(round(slope,4)) + ' \n' +' \n' +'NPQmax = ' + str(round(NPQmax,4))
                          +' \n' +  'Offest(x,y) =' + str(round(PARsv_offset_X18_S08,2)) + ',' + str(round(NPQsv_offset_X18_S08,2) ) 
                          , xy=(0, 1), xytext=(12, -12), va='top', xycoords='axes fraction', textcoords='offset points', color='k')
            
            plt.legend(loc='upper right')


            plt.suptitle ('Profil vertical fluo {0}'.format(data['CYCLE_NUMBER'].values.astype(int)[i]) +' {0} '.format(dire)+WMO + '  '+ d)



# =============================================================================
#%% Other Fig

# Visualize NPQ versus depth
            # plt.figure()
            # plt.plot(NPQsv_plot,depth_NPQ,label='NPQsv',color='grey')
            # plt.xlabel('NPQsv')
            # plt.ylabel('Depth')
            # plt.gca().invert_yaxis()
            
            # plt.twiny()
            # plt.plot(PARsv_plot,depth_NPQ,label='iPARsv',color='red')
            
            # plt.title('Vision de NPQ face au PAR')
            
# Temporal serie of alpha NPQ and PAR mean in MLD

            # plt.figure(figsize=[12,8])
           
            
            # plt.subplot(211)
            # plt.title('Evolution temporelle de '+ r'$\alpha$_NPQ' + ' flotteur WMO:' + WMO)
            # plt.plot(time_NPQ_X18_S08,alpha_NPQ_X18_S08,'b',label='Alpha_NPQ')
            # plt.ylabel(r'$\alpha$_NPQ')
            # plt.xticks(time_NPQ_X18_S08[::int(len(time_NPQ_X18_S08)/10)])
            # plt.legend()
           
            # from datetime import datetime, timedelta
            # from scipy.interpolate import interp1d
            
            
            
            # # fit_alpha = interp1d(time_NPQ_X18_S08, alpha_NPQ_X18_S08, kind='cubic')
            # # t = np.arange(datetime(2015,5,10), datetime(2018,11,18), timedelta(days=1)).astype(datetime)
            # # Soca_Res05 = fit_alpha(t)
            
            # # plt.plot(t,Soca_Res05)
            
            
            # # plt.twinx()
            # plt.subplot(212)
            # plt.plot(time_NPQ_X18_S08,PAR_mean_X18_S08,'k--',label='MLD')
            # plt.ylabel('Depth')
            # plt.xticks(time_NPQ_X18_S08[::int(len(time_NPQ_X18_S08)/10)])
            # plt.legend()
            # # plt.gca().invert_yaxis()
#%% Save array

    if  len(alpha_NPQ_X18_S08) != 0: # verifie si au moins 1 profile valide
        
            
            print(str(WMO)+': ' + str(len(alpha_NPQ_X18_S08)))
           
           
            
            alpha_NPQ_txt_X18_S08.extend(alpha_NPQ_X18_S08)
            
            
            time_NPQ_txt_X18_S08.extend(time_NPQ_X18_S08) 
            
            latitude_txt_X18_S08.extend(latitude_X18_S08)
            longitude_txt_X18_S08.extend(longitude_X18_S08)
            
            MLD_txt_X18_S08.extend(MLD_NPQ_X18_S08)
            PAR_mean_txt_X18_S08.extend(PAR_mean_X18_S08)



# %% Ecriture du .txt (csv max row = 1000)

# data = {"Time" : time_NPQ_txt_X18_S08 , "Longitude" : longitude_txt_X18_S08 , "Latitude":latitude_txt_X18_S08 ,
#         "alpha" :alpha_NPQ_txt_X18_S08 ,"MLD" : MLD_txt_X18_S08 , "PAR_moyen" : PAR_mean_txt_X18_S08} 
# data = pd.DataFrame(data) 

# data.to_csv('/home/oao2020/Documents/Pierrick/Outputs/alpha_NPQ_X18_S08_ALL_CHL_DO_V5.txt', header=True, index=None, sep=';', mode='a')







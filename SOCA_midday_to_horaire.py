# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 19:43:44 2024

Objectif : Création d'un NetCDF par flotteur où les données des profils climatologiques et synhtétiques de PAR simulés par SOCA sont stockées si disponibles

Entrée : Profils de PAR simulé par SOCA au midi local
- Construction d'un profil heure locale : PAR de surface produit en croix depuis une cloche d'éclairement (pvlib) 
                                          et KdPAR(z) depuis le profil initial
Sortie : NetCDF


Rejet de profil : besoin que mesure CHLa, et que QC physique soit bons


@author: OMTABmezz
"""
# =============================================================================
import warnings
# Ignore runtime warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
# =============================================================================

import os

import xarray as xr
import pandas as pd
import numpy as np



import timezonefinder, pytz
from datetime import datetime

import gsw







#%% Functions

def MLD_Func (x,y):
    """
    x : Densité
    y : Profondeur
    
    seuil : densité à 10m + 0.03kg/m3 (BM 2004)

    """
    y = y[~np.isnan(x)] # garde que y associé a une valeur de x
    x = x[~np.isnan(x)]
   
    index = np.absolute(y - 10).argmin() #mieux prend profondeur plus proche
    # index = np.max(np.where(y < 10)) +1 # index d'environ 10dBar
    dix = x[index] #valeur à 10m
    
    
    index = np.min(np.where(x[index:] > dix + 0.03 ) +index) #faire plus index car faire [:] raccourci la liste
    MLD = y[index]
    
    return MLD

def time_UTCtolocal (time_64,longitude):
    
    delta = 60*longitude /15  # 60 min = 15°
    
    local = time_64 + np.timedelta64(int(delta),'m')
    
    return local


from pvlib.location import Location
# https://pvlib-python.readthedocs.io/en/stable/user_guide/clearsky.html

def Radiance_from_midday (date, LAT, LON, Radiance_LocalNoon, fuseau_horaire) :
    """
    
    Parameters
    ----------
    date : pd.datetime
        Local time 
        
    LAT : float
    LON : float
    
    Radiance_LocalNoon : float
        Solar irrandiance at midday [W/m2] (or PAR)

    Returns Solar irradiance (or PAR) at an other hour x
    -------
     /!\ tz dans position est par défault UTC et donc fait des décalage car je travail en local moi 
    """
    # Position info
    position = Location(latitude=LAT, longitude=LON, altitude = 0, tz = fuseau_horaire)
    
    # Time : pandas.DatetimeIndex, LOCAL time !
    year = date.year
    month = date.month
    day = date.day
    heure_local = date.hour
    
    start_day = str(year) + '-' + str(month) + '-' + str(day)+ ' 00:00:00'
    end_day = str(year) + '-' + str(month) + '-' + str(day)+ ' 23:30:00'
    
    # Time array 
    time = pd.date_range(start = start_day, end = end_day, freq='30min', tz = fuseau_horaire)
    
    # Gaussian cruve of GHI = Global horizontal irrandiance
    radiance = position.get_clearsky(times = time, model = 'ineichen')
    
    
    # Get solar irrandiance at another hour
    # GHI = Global horizontal irradiance = Direct normal irradiance (dni) + Diffuse horizontal irradiance (dhi) + cos(z)
    GHI_midday = radiance.loc[str(year) + '-' + str(month) + '-' + str(day)+ ' 12:00:00', 'ghi']
    
    index_nearest = np.argmin(np.abs(time - date.tz_localize(fuseau_horaire)))
    GHI_xH = radiance.loc[time[index_nearest], 'ghi']
    # GHI_xH = radiance.loc[str(year) + '-' + str(month) + '-' + str(day)+ ' '+str(heure_local)+':00:00', 'ghi']
    
    Radiance_xH = Radiance_LocalNoon * GHI_xH / GHI_midday
    
    
    return Radiance_xH, radiance




def PAR_zt_frm_PAR_midday (date, LAT, LON, PAR_SOCAcsv_LocalNoon) :
    """
    Bien donner date en local ! et sous forme pd to datetime

    """
    try :
        # if dataframe
        PAR_midi = PAR_SOCAcsv_LocalNoon.PAR_ADJUSTED[0]
        
        tf = timezonefinder.TimezoneFinder()
        fuseau_horaire = tf.timezone_at(lat=LAT, lng=LON)
        
        PAR_xH , radiance = Radiance_from_midday (date , LAT, LON, PAR_midi, fuseau_horaire)
        
            
        build_PAR = np.array([])
        build_PAR = np.append(build_PAR, PAR_xH)
        
        for i in range(len(PAR_SOCAcsv_LocalNoon.PAR_ADJUSTED.values)-1) : 
            
            PAR_z = build_PAR[i] * np.exp(-np.abs(np.diff(np.round(np.log(PAR_SOCAcsv_LocalNoon.PAR_ADJUSTED.dropna().values),5))/5)[i] * 5) # 5 = en mètres = résolution de SOCA
            build_PAR = np.append(build_PAR, PAR_z)
            
            
    except :
        # if array
        PAR_midi = PAR_SOCAcsv_LocalNoon[0]
    
        tf = timezonefinder.TimezoneFinder()
        fuseau_horaire = tf.timezone_at(lat=LAT, lng=LON)
        
        PAR_xH , radiance = Radiance_from_midday (date , LAT, LON, PAR_midi, fuseau_horaire)
        
            
        build_PAR = np.array([])
        build_PAR = np.append(build_PAR, PAR_xH)
    
        for i in range(len(PAR_SOCAcsv_LocalNoon)-1) : 
            
            PAR_z = build_PAR[i] * np.exp(-np.abs(np.diff(np.round(np.log(PAR_SOCAcsv_LocalNoon),5))/5)[i] * 5) # 5 = en mètres = résolution de SOCA
            build_PAR = np.append(build_PAR, PAR_z)
    
    return build_PAR


    
#%% Loop
 
# Read Sprof_files path , floats with CHLA
df = pd.read_csv("D:/DATA/data_flotteurs_CSV/Stage M2/Sprof_sorted_01_2024.csv")
# Le faire pour tous les flotteurs mesurant de la CHL
column = 'CHL'


# =============================================================================
path_climato = 'D:/DATA/data_flotteurs_CSV/Stage M2/SOCA_light/SOCA_LIGHT_Weekly_Climatology/SOCA_LIGHT_Weekly_Climatology/20240115/'
# Chargement en amont des NETCDF climatologiques (normalement devrait être plus rapide)
Climato_dataset = { }
for semaine in np.arange(1,52+1,1) :
        numero_semaine = str(semaine).zfill(2) 
        
        fichier = os.listdir(path_climato + numero_semaine)[0]
        Climato_dataset[numero_semaine] = xr.open_dataset(path_climato + str(numero_semaine) +'/' + fichier)
# =============================================================================
        

# Fonction pour tentative de parralélisation avec concurrent future
def Save_SOCA_NetCDF (Sprof_path) :
    
    # Ouverture du fichier Sprof
    try : 
        data = xr.open_dataset(Sprof_path)
    except : 
        # Path to CHLA_no BBP
        new_path = '/'.join(Sprof_path.split('/')[:3]) + '/CHLA_but_noBBP/' + '/'.join(Sprof_path.split('/')[3:])
        data = xr.open_dataset(new_path)
# =============================================================================   
    # Extract some variables
    time = data['JULD'].values 
    LON = data['LONGITUDE'].values 
    LAT = data['LATITUDE'].values 
    
    WMO = Sprof_path.split('/')[-1].split('_')[0]
    # dac = Sprof_path.split('/')[-2]
    CYCLE_NUMBER = data['CYCLE_NUMBER'].values.astype(int)
    
# =============================================================================
#   Création du NCDF remplit de NAN values
# =============================================================================
    # Dimensions
    cycle_number_dim = len(CYCLE_NUMBER)
    depth_dim = len(np.arange(0, 255, 5))
    
    ds = xr.Dataset({
    'PAR_SOCA': ([ 'CYCLE_NUMBER', 'depth'], np.full((cycle_number_dim , depth_dim), np.nan)),
    'PAR_climato': (['CYCLE_NUMBER', 'depth'], np.full(( cycle_number_dim, depth_dim ), np.nan)),},
        
    coords={ 'CYCLE_NUMBER':  CYCLE_NUMBER,
            'depth':  np.arange(0, 255, 5)})
# =============================================================================    
    # Output path
    NC_Output_path = 'D:/DATA/data_flotteurs_CSV/Stage M2/SOCA_light/PAR_Outputs/Midday_to_Horaire_py/' + str(WMO) + '_PAR_SOCA_climato_cloche.nc'
    print(WMO)
    # Check si fichier existe
    if not os.path.exists(NC_Output_path) :
      
# =============================================================================
# Loop on cycle number
# =============================================================================

        for i in range(len(data['CYCLE_NUMBER'].values)): 
            
            # Pour ne pas avoir à faire de continue = problèmes avec threadpool
            # continue_bool = False
            
            # Déterminer l'heure locale
            try :
                ## Méthode 1 : avec module timezone = plus gourmand en calcul
                # tf = timezonefinder.TimezoneFinder()
                # timezone = pytz.timezone(tf.certain_timezone_at(lat = LAT[i], lng = LON[i]))
                # date_local = pd.to_datetime(time[i]) + timezone.utcoffset(pd.to_datetime(time[i]))
                
                ## Méthode 2 : Ma fonction UTC to local avec numpy : bien plus rapide et résultats proches 
                date_local =  pd.to_datetime(time_UTCtolocal(time[i],LON[i]))
                
                # Avoir le format de la date, permet aussi de détecter les nan si time[i] = Nat
                date_local = pd.to_datetime(date_local.strftime('%Y-%m-%d %H:%M:%S'))
            except :
                print('Long or time NaN ' + str(i))
                # continue_bool = True 
                continue
            

            # if not date_local.hour in np.arange(7,18+1,1) :
            #     print('nuit')
            
            # test si dans tranche horaire 7-18h  
            if date_local.hour in np.arange(7,18+1,1):   
                
                print(str(int(i)) + '/' + str(int(np.unique(data['CYCLE_NUMBER'].values)[-1])))


#%%  PHYSIQUE Conditions
# =============================================================================
# Prise Adjusted si présent et si plus de 75% de QC = 1,2,5,8 (Si QC mauavais pas test du RT car sera forcément pas bon)
# Si Adjusted Absent : 
# Prise de RT et check si bon QC > 75%      
# =============================================================================
                depth = data['PRES_ADJUSTED'].values[i]
                
                if not any(~np.isnan(depth)):
                    
                    depth = data['PRES'].values[i]
                    
                    if not any(~np.isnan(depth)):
                        print(str(i) +'Pas de mesure de Pression Rt ni ADJ ' )
                        continue
                        # continue_bool = True
                    else:
                    
                         depth_QC = data['PRES_QC'].values[i]
                         depth_QC = depth_QC[~pd.isnull(depth_QC)].astype('int')  
                         profile_QC = np.where ((depth_QC ==1) | (depth_QC==2)|(depth_QC ==5) | (depth_QC==8))
                         percent_depth = 100* len(profile_QC[0])/len(depth_QC)
                         
                         if not percent_depth > 75:
                           print(str(i) +': QC profil PRES trop mauvais : '+ str(percent_depth))
                           continue
                           # continue_bool = True
                else : 
                     depth_QC = data['PRES_ADJUSTED_QC'].values[i]
                     depth_QC = depth_QC[~pd.isnull(depth_QC)].astype('int')  
                     profile_QC = np.where ((depth_QC ==1) | (depth_QC==2)|(depth_QC ==5) | (depth_QC==8))
                     percent_depth = 100* len(profile_QC[0])/len(depth_QC)
                     
                     if not percent_depth > 75:
                       print(str(i) +': QC profil PRES_ADJ trop mauvais : '+ str(percent_depth))
                       continue
                       # continue_bool = True
          
        # =============================================================================
         
                SP = data['PSAL_ADJUSTED'].values[i]
                
                if not any(~np.isnan(SP)):
                    
                    SP = data['PSAL'].values[i]
                    
                    if not any(~np.isnan(SP)):
                        print(str(i) +'Pas de mesure de PSAL Rt ni ADJ ' )
                        continue
                        # continue_bool = True
                    
                    else:
                        SP_QC = data['PSAL_QC'].values[i]
                        SP_QC = SP_QC[~pd.isnull(SP_QC)].astype('int')  
                        profile_QC = np.where ((SP_QC ==1) | (SP_QC==2)|(SP_QC ==5) | (SP_QC==8))
                        percent_SP = 100* len(profile_QC[0])/len(SP_QC)
                        
                        if not percent_SP > 75:
                          print(str(i) +': QC profil SP trop mauvais : '+ str(percent_SP))
                          continue
                          # continue_bool = True
                        
                else : 
                     SP_QC = data['PSAL_ADJUSTED_QC'].values[i]
                     SP_QC = SP_QC[~pd.isnull(SP_QC)].astype('int')  
                     profile_QC = np.where ((SP_QC ==1) | (SP_QC==2)|(SP_QC ==5) | (SP_QC==8))
                     percent_SP = 100* len(profile_QC[0])/len(SP_QC)
                     
                     if not percent_SP > 75:
                       print(str(i) +': QC profil SP ADJ trop mauvais : '+ str(percent_SP))
                       continue
                       # continue_bool = True
                
                      
        # =============================================================================      
                
                temp = data['TEMP_ADJUSTED'].values[i]
                
                if not any(~np.isnan(temp)):
                    
                    temp = data['TEMP'].values[i]
                    
                    if not any(~np.isnan(temp)):
                        print(str(i) +'Pas de mesure de TEMP Rt ni ADJ ' )
                        continue
                        # continue_bool = True
                    else:
                        Temp_QC = data['TEMP_QC'].values[i]
                        Temp_QC = Temp_QC[~pd.isnull(Temp_QC)].astype('int')  
                        profile_QC = np.where ((Temp_QC ==1) | (Temp_QC==2)|(Temp_QC ==5) | (Temp_QC==8))
                        percent_Temp = 100* len(profile_QC[0])/len(Temp_QC)
                        
                        if not percent_Temp > 75:
                          print(str(i) +': QC profil TEMP trop mauvais : '+ str(percent_Temp))
                          continue
                          # continue_bool = True
                        
                else : 
                     Temp_QC = data['TEMP_ADJUSTED_QC'].values[i]
                     Temp_QC = Temp_QC[~pd.isnull(Temp_QC)].astype('int')  
                     profile_QC = np.where ((Temp_QC ==1) | (Temp_QC==2)|(Temp_QC ==5) | (Temp_QC==8))
                     percent_Temp = 100* len(profile_QC[0])/len(Temp_QC)
                     
                     if not percent_Temp > 75:
                       print(str(i) +': QC profil TEMP trop mauvais : '+ str(percent_Temp))
                       continue
                       # continue_bool = True
                
        
    # =============================================================================
    #    Compute density and find MLD
    # =============================================================================
            
                SA = gsw.SA_from_SP(SP,depth,LON[i],LAT[i])
                CT = gsw.CT_from_t(SA,temp,depth)
                # Calcul de la densité avec le protocole TEOS10 (salinité absolue et température conservative)
                density = gsw.sigma0(SA,CT)
                
                # Rm Nan values
                depth_dens = depth[~np.isnan(density)]
                density = density[~np.isnan(density)]
    
    
                if not any(depth_dens<10) or all(depth_dens<10):
                    print(str(i) +' Pas de physique en surface : depth_dens[0] = ' + str(depth_dens[:0]))
                    continue
                    # continue_bool = True
                
    
                try :
                    MLD = MLD_Func(density, depth_dens)
                except :
                    print('Soucis MLD')
                    continue
                    # continue_bool = True
                
                
                # Si l'entrée n'a pas été uptade ==> dans le cas de parralelisation ou n'accepte pas continue
                # if continue_bool == False : 


            
#%%  Lecture PAR SOCA           
# =============================================================================
# Lecture PAR SOCA
# ============================================================================= 
                # Lecture des fichiers à midi local
                path_SOCA_outputs = 'D:/DATA/data_flotteurs_CSV/Stage M2/SOCA_light/PAR_Outputs/Midday/' +  str(WMO) + '/'
        
                direction = data['DIRECTION'].astype(str).values[i]
                try : 
                    if direction == 'A' :
                        SOCA = pd.read_csv(path_SOCA_outputs  + WMO + '_' + str(CYCLE_NUMBER[i]) + '_PAR_12h.txt' )
                    else :
                        SOCA = pd.read_csv(path_SOCA_outputs  + WMO + '_' + str(CYCLE_NUMBER[i]) + 'D_PAR_12h.txt' )
                                                   
# =============================================================================
# Reconstruire un profil SOCA depuis une valeur de surface à midi et les KdPAR 
# =============================================================================
                                                    
                    PAR_SOCA =  PAR_zt_frm_PAR_midday (date_local, LAT[i], LON[i], SOCA)    
                    
                except:
                    print('Pas de fichier PAR_SOCA pour ce profil : '+WMO + '_' + str(CYCLE_NUMBER[i]) + direction + '_PAR')
                    PAR_SOCA = np.full(51, np.nan)
        

                # Save dans le Nc
                ds['PAR_SOCA'][i,:] = PAR_SOCA

#%% Lecture PAR Climato
# =============================================================================
# Lecture PAR Climato
# =============================================================================
        
                # Ouverture du fichier
                # 1/ Recherche du numero de la semaine
                date_object = datetime.strptime(str(date_local), '%Y-%m-%d %H:%M:%S')
                numero_semaine = date_object.strftime('%W')
                # Si année bissextile donne 00 le 1er janvier et 53 le 31 décembre
                if numero_semaine == '00' : 
                    numero_semaine = '01'
                    
                elif numero_semaine == '53' :
                    numero_semaine = '52'
                
                # Prise du Nc préalablement stocké en mémoire vive
                climato = Climato_dataset[numero_semaine]
    
                
                # recherche de la position la plus proche de celle du flotteur
                index_lon = np.absolute(LON[i] - climato.longitude.values ).argmin()
                index_lat = np.absolute(LAT[i] - climato.latitude.values ).argmin()
                
                PAR_climato = climato.PAR[0,:,index_lat,index_lon].values
    
# =============================================================================
# La climato hebdomadaire n'est pas remplit à 100% comme l'était celle mensuelle
# Si pas de données : recherche 2 longitudes les plus proches (hypothèse moins de variations longitunidalement)
# Si toujours pas recherche sur les latitudes les plus proches 

                if np.all(np.isnan(PAR_climato)) :
                    # Recherche des 3 longitudes les plus proches de la position flotteur
                    Long_nearest = np.argsort(np.absolute(LON[i] - climato.longitude.values ))[:3]
                    # Recherche du PAR a ces localisation
                    near_1 = climato.PAR[0,:,index_lat , Long_nearest[1]].values
                    near_2 = climato.PAR[0,:,index_lat , Long_nearest[2]].values
                    # Faire une moyenne des deux 
                    PAR_climato = np.nanmean(np.vstack((near_1, near_2)), axis = 0)
                    
                    if np.all(np.isnan(PAR_climato)) :
                        # Si aucune donnée sur la même longitude test de la même latitude
                        # Recherche des 3 latitudes les plus proches de la position flotteur
                        Lat_nearest = np.argsort(np.absolute(LAT[i] - climato.latitude.values ))[:3]
                        near_1 = climato.PAR[0,:,Lat_nearest[1] ,index_lon ].values
                        near_2 = climato.PAR[0,:,Lat_nearest[2] , index_lon ].values
                        
                        PAR_climato = np.nanmean(np.vstack((near_1, near_2)), axis = 0)
                        
                        if np.all(np.isnan(PAR_climato)) :
                            print('Pas de données climato')
                            # continue # Pas besoin car de toute manière sont deja remplit de nan
                    
# =============================================================================
# =============================================================================
# Reconstruire un profil SOCA depuis une valeur de surface à midi et les KdPAR
# Et save to csv 
# =============================================================================

                PAR_climato =  PAR_zt_frm_PAR_midday (date_local, LAT[i], LON[i], PAR_climato)
        
                # Save dans le Nc
                ds['PAR_climato'][i,:] = PAR_climato
      
                
        return ds.to_netcdf(NC_Output_path)



#%% Boucle FOR


# for Sprof_path in df['CHL'] : 
#     Save_SOCA_NetCDF(Sprof_path)



#%% Parallèlisation : Marche pas de fou jsp pourquoi encore 

# Hypothèse : python bloque la parralilation si plusieurs ouvriers stocks des données aux mêmes endroits
# Et donc potentiellement plus lent qu'une simple boucle for
# ==> faire des recherche sur ProcessPoolExecutor, pas paraléllisation en Thread mais par processus 

# import concurrent.futures

# def download_files_in_parallel(Sprof_path):

            
        
#         with concurrent.futures.ThreadPoolExecutor() as executor:
#             futures = []
            
            
            
#             for row in Sprof_path: 
                
#                 future = executor.submit(Save_SOCA_NetCDF, row)
#                 futures.append(future)
                  
            
#             # Attendre la fin de tous les téléchargements
#             concurrent.futures.wait(futures)


# # =============================================================================
# # =============================================================================

# Sprof_path = df['CHL']

# download_files_in_parallel(Sprof_path[:3])








   
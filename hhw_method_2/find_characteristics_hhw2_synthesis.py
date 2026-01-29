#!/usr/bin/env python
# coding: utf-8


# Create a dataframe of the characteristics of observed HHW2


# In[1]:


import numpy as np
import pandas as pd


# In[2]:


dirin       = '/OBS/HIM/'
dirres= "/OBS/heatwaves/"


# In[3]:


Stations=["Gillot","Pamandzi","Tromelin"]


# In[4]:


for station in Stations:
        
        
        data_hist = np.load(dirin+station+'_heatwaves_HIM_obs.npy')
        
        df=pd.DataFrame(columns=["year","length","intensity","severity"])
        
                
        
        df["year"]=data_hist[:,0].astype(int)
        df["length"]=data_hist[:,1].astype(int)
        df["intensity"]=data_hist[:,2]
        df["severity"]=data_hist[:,3]#severity is not described in the paper
        
        df.to_csv(dirres+"/"+station+"_carac_heatwaves_DATAFRAME_HIM_1990-2014.csv")        
        
        


################################################################################

# Create a dataframe of the characteristics of past and future HHW2 in CMIP6 models


# In[1]:


import numpy as np
import pandas as pd


# In[2]:


dirin = '/CMIP6/HIM/'
dirres= "/CMIP6/heatwaves/"


# In[3]:


Modeles=["ACCESS-CM2",
         "CanESM5",
         "CNRM-ESM2-1",
    "CNRM-CM6-1",
    "CNRM-CM6-1-HR",
         "EC-Earth3","EC-Earth3-CC","EC-Earth3-Veg-LR",
         "KACE-1-0-G",
         "IPSL-CM6A-LR","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR",
         "MPI-ESM1-2-LR",
         "NorESM2-MM"
        ]


Stations=["Gillot","Pamandzi","Tromelin"]


# In[4]:


for station in Stations:
    
    print(station+'\n')
    
    df_station_ssp245=pd.DataFrame(columns=["year","length","intensity","severity"])
    df_station_ssp585=pd.DataFrame(columns=["year","length","intensity","severity"])
    # severity is not used in the paper
    
    for modele in Modeles:
        
        print(modele)
        
        hw_hist = np.load(dirin+station+'_'+modele+'_heatwaves_HIM_hist.npy')        
        hw_ssp245 = np.load(dirin+station+'_'+modele+'_heatwaves_HIM_fut_ssp245.npy')
        hw_ssp585 = np.load(dirin+station+'_'+modele+'_heatwaves_HIM_fut_ssp585.npy')
        
        
        if hw_hist.shape[0]==0:
            
            data_ssp245=hw_ssp245
            data_ssp585=hw_ssp585
        else:
                    
            data_ssp245=np.concatenate((hw_hist,hw_ssp245),axis=0)
            data_ssp585=np.concatenate((hw_hist,hw_ssp585),axis=0)
        
        df_ssp245=pd.DataFrame()
        df_ssp585=pd.DataFrame()
        
        df_ssp245["year"]=data_ssp245[:,0].astype(int)
        df_ssp245["length"]=data_ssp245[:,1].astype(int)
        df_ssp245["intensity"]=data_ssp245[:,2]
        df_ssp245["severity"]=data_ssp245[:,3]
        
        df_ssp585["year"]=data_ssp585[:,0].astype(int)
        df_ssp585["length"]=data_ssp585[:,1].astype(int)
        df_ssp585["intensity"]=data_ssp585[:,2]
        df_ssp585["severity"]=data_ssp585[:,3]
                
        df_ssp245.to_csv(dirres+"/"+station+"_"+modele+"_carac_heatwaves_DATAFRAME_HIM_ssp245_1990-2100.csv")
        df_ssp585.to_csv(dirres+"/"+station+"_"+modele+"_carac_heatwaves_DATAFRAME_HIM_ssp585_1990-2100.csv")
        
        df_station_ssp245=pd.concat([df_station_ssp245,df_ssp245])
        df_station_ssp585=pd.concat([df_station_ssp585,df_ssp585])
        
            
    df_station_ssp245.to_csv(dirres+"/"+station+"_ALL_MODELES_carac_heatwaves_DATAFRAME_HIM_ssp245_1990-2100.csv")
    df_station_ssp585.to_csv(dirres+"/"+station+"_ALL_MODELES_carac_heatwaves_DATAFRAME_HIM_ssp585_1990-2100.csv") 




################################################################################


# Create a synthesis of CMIP6 and ALADIN HHW2


# In[1]:


import numpy as np
import pandas as pd


# In[2]:


dirres= "/CMIP6/heatwaves/"


# In[3]:


Modeles=["ACCESS-CM2",
         "CanESM5",
         "CNRM-ESM2-1",
    "CNRM-CM6-1",
    "CNRM-CM6-1-HR",
         "EC-Earth3","EC-Earth3-CC","EC-Earth3-Veg-LR",
         "KACE-1-0-G",
         "IPSL-CM6A-LR","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR",
         "MPI-ESM1-2-LR",
         "NorESM2-MM"
        ]

range_date=np.arange(1995,2096)

Stations=["Gillot","Pamandzi","Tromelin"]

Scenarios=["ssp245","ssp585"]


# In[4]:


def calcul_carac_heatwaves(station,scenario):
    
    '''
    This method create files that summarizes the HHW2 characteristics (duration, intensity and severity)
    for CMIP6 and ALADIN.
    We get results sorted by decades, these decades covering the entire period 1990-2100.
    '''
    
    
    i=0
    
    #CMIP6
        
    data_station=pd.read_csv(dirres+"/"+station+"_ALL_MODELES_carac_heatwaves_DATAFRAME_HIM_"+scenario+"_1990-2100.csv",header=0,names=['year','length','intensity','severity'])
    
    data_station_sort=data_station.sort_values(by='year')
    
    
    if station=="Gillot":
    
        #RCM ALADIN -- Only for Gillot
    
        dirres_aladin="/ALADIN/heatwaves/"
    
        data_station_aladin=pd.read_csv(dirres_aladin+station+"_carac_heatwaves_DATAFRAME_HIM_"+scenario+"_1990-2100.csv",header=0,names=['year','length','intensity','severity'])
    
        data_station_sort_aladin=data_station_aladin.sort_values(by='year')

    
    dec=[]
    
    for year in data_station_sort.year.values:
        
        if not year == 1990:
            
            dec_year=((year-1)//10)*10 + 5
        
            dec.append(dec_year)
            
        else:
            
            dec.append(1995)
            
    data_station_sort["dec"]=dec

    
    if station == "Gillot":
        dec=[]
    
        for year in data_station_sort_aladin.year.values:
        
            if not year == 1990:
            
                dec_year=((year-1)//10)*10 + 5
        
                dec.append(dec_year)
            
            else:
            
                dec.append(1995)
            
    
    
        data_station_sort_aladin["dec"]=dec
    
    
    for dec_value in np.unique(data_station_sort.dec.values):
        
        
        data_dec=data_station_sort[data_station_sort.dec==dec_value]
        
        np.save(dirres+station+"_length_heatwaves_"+scenario+"_CMIP6_"+str(dec_value)+".npy",data_dec.length.values)   
        np.save(dirres+station+"_intensity_heatwaves_"+scenario+"_CMIP6_"+str(dec_value)+".npy",data_dec.intensity.values)    
        np.save(dirres+station+"_severity_heatwaves_"+scenario+"_CMIP6_"+str(dec_value)+".npy",data_dec.severity.values) 
        
        
        if station=='Gillot': #add RCM ALADIN for Gillot
            data_dec_aladin=data_station_sort_aladin[data_station_sort_aladin.dec==dec_value]
        
            np.save(dirres_aladin+station+"_length_heatwaves_"+scenario+"_ALADIN_"+str(dec_value)+".npy",data_dec_aladin.length.values)   
            np.save(dirres_aladin+station+"_intensity_heatwaves_"+scenario+"_ALADIN_"+str(dec_value)+".npy",data_dec_aladin.intensity.values)    
            np.save(dirres_aladin+station+"_severity_heatwaves_"+scenario+"_ALADIN_"+str(dec_value)+".npy",data_dec_aladin.severity.values)    


# In[5]:


for station in Stations:
    for scenario in Scenarios:
        calcul_carac_heatwaves(station,scenario)


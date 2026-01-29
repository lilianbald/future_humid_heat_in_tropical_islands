#!/usr/bin/env python
# coding: utf-8

# ## Content: Compute CMIP6 uncertainties for number of days/year under HHW2 conditions

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


dirres= "/CMIP6/heatwaves/"
dirobs="/OBS/heatwaves/"


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


def bootstrap(X,n):
    res=np.empty((n,X.shape[1]))
    randind=np.random.choice(np.arange(X.shape[0]),(n,))
    for i in range(n):
        res[i,:]=X[randind[i],:]
    return res


# In[5]:


def calcul_nb_days_heatwaves(station,scenario):
    
    '''
    Return quantile 0.1, 0.9, median and all bootstrap samples
    Data were formerly smoothed over a 10-year sliding window
    '''
    
    data_stations=np.empty((len(Modeles),len(range_date)))
    data_stations_hw=np.empty((len(Modeles),len(range_date)))

    
    i=0
    
    
    for modele in Modeles:
        
        data_station=pd.read_csv(dirres+"/"+station+"_"+modele+"_smoothed_heatwave_days_"+scenario+".csv",header=0,names=['year','nb_days'])
        
        data_stations[i,:]=data_station.nb_days.values[:]
        
        data_station_hw=pd.read_csv(dirres+"/"+station+"_"+modele+"_smoothed_heatwave_days_over3days_"+scenario+".csv",header=0,names=['year','nb_days'])
        
        data_stations_hw[i,:]=data_station_hw.nb_days.values[:]
        
        i+=1
        
    n=1000000
    
    
    samples=bootstrap(data_stations,n)
    
    q10 = np.nanquantile(samples,0.10,axis=0)
    q90 = np.nanquantile(samples,0.90,axis=0)
    median= np.nanmedian(samples,axis=0)
    
    q10_df=pd.DataFrame(index=range_date.tolist())
    q10_df['nb_days']=q10
    q90_df=pd.DataFrame(index=range_date.tolist())
    q90_df['nb_days']=q90
    median_df=pd.DataFrame(index=range_date.tolist())
    median_df['nb_days']=median
    
    q10_df.to_csv(dirres+station+"_q10_smoothed_heatwave_HIM_days_"+scenario+"_CMIP6.csv")
    q90_df.to_csv(dirres+station+"_q90_smoothed_heatwave_HIM_days_"+scenario+"_CMIP6.csv")
    median_df.to_csv(dirres+station+"_q50_smoothed_heatwave_HIM_days_"+scenario+"_CMIP6.csv")
    np.save(dirres+station+"_samples_smoothed_heatwave_HIM_days_"+scenario+"_CMIP6.npy",samples)
    
    
    samples=bootstrap(data_stations_hw,n)
    
    q10 = np.nanquantile(samples,0.10,axis=0)
    q90 = np.nanquantile(samples,0.90,axis=0)
    median= np.nanmedian(samples,axis=0)
    
    q10_df=pd.DataFrame(index=range_date.tolist())
    q10_df['nb_days']=q10
    q90_df=pd.DataFrame(index=range_date.tolist())
    q90_df['nb_days']=q90
    median_df=pd.DataFrame(index=range_date.tolist())
    median_df['nb_days']=median
    
    q10_df.to_csv(dirres++station+"_q10_smoothed_heatwave_HIM_days_over3days_"+scenario+"_CMIP6.csv")
    q90_df.to_csv(dirres+station+"_q90_smoothed_heatwave_HIM_days_over3days_"+scenario+"_CMIP6.csv")
    median_df.to_csv(dirres++station+"_q50_smoothed_heatwave_HIM_days_over3days_"+scenario+"_CMIP6.csv")
    np.save(dirres+station+"_samples_smoothed_heatwave_HIM_days_over3days_"+scenario+"_CMIP6.npy",samples)
    
    


# In[6]:


for station in Stations:
    for scenario in Scenarios:
           calcul_nb_days_heatwaves(station,scenario)


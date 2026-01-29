#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


Stations=["Gillot","Pamandzi","Tromelin"]


# In[3]:


def count_periods_over3days(data,thres):
    
    '''
    This method aims at counting days under HHW1 conditions in binary observation series  (1 if HIX>=s, 0 otherwise)
    thresholds:
    s1 : CAUTION
    s2 : DANGER
    s3 : EXTREME DANGER
    '''
    
    s=data[thres].values.copy()

    times=data.time.values
    
    i=0
    
    j=0

    
    while i<(times.shape[0]):
        
        
        if s[i]==1:
            
            j+=1
            
        
        else:
            
            if j<3 and j>0 and i>2:
                
                s[i-1]=0
                s[i-2]=0
            
            j=0
        
            
        i+=1
    
    df=pd.DataFrame(columns=["time",thres])
    
    df['time']=pd.to_datetime(data.time)
    df[thres]=s
    
    data['time']=pd.to_datetime(data.time)
    
    data_year=data[thres].groupby(data.time.dt.year).sum()
    df_year=df[thres].groupby(df.time.dt.year).sum()
    
    
    df.to_csv("/OBS/heatwaves/"+station+"_OBS_binary_NOAA_"+thres+"_HIX_over3days_1985-2014.csv")
    df_year.to_csv("/OBS/heatwaves/"+station+"_OBS_nb_days_over3days_NOAA_"+thres+"_HIX_1985-2014.csv")

    
    
    


# In[4]:


for station in Stations:

    # file with: 1 if HIX>=thres, 0 otherwise
    obs_binary=pd.read_csv("/OBS/heatwaves/"+station+"_OBS_binary_NOAA_thresholds_HIX_1985-2014.csv")
    
    for thres in ['s1','s2','s3']:
        
        count_periods_over3days(obs_binary,thres)
        
    

#######################################################################################


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


Stations=["Gillot","Pamandzi","Tromelin"]


# In[4]:


def count_periods_over3days(data_ssp245,data_ssp585,thres):
    
    '''
    This method aims at counting days under HHW1 conditions in binary CMIP6 series  (1 if HIX>=s, 0 otherwise)
    thresholds:
    s1 : CAUTION
    s2 : DANGER
    s3 : EXTREME DANGER
    Values are smoothed over 10 years
    '''
    
    s_ssp245=data_ssp245[thres].values.copy()
    s_ssp585=data_ssp585[thres].values.copy()

    times_ssp245=data_ssp245.time.values
    
    i=0
    
    j=0

    
    while i<(times_ssp245.shape[0]):
        
        
        if s_ssp245[i]==1:
            
            j+=1
            
        
        else:
            
            if j<3 and j>0 and i>2:
                
                s_ssp245[i-1]=0
                s_ssp245[i-2]=0
            
            j=0
        
            
        i+=1
    
    df_ssp245=pd.DataFrame(columns=["time",thres])
    
    df_ssp245['time']=pd.to_datetime(data_ssp245.time)
    df_ssp245[thres]=s_ssp245

    df_year_ssp245=df_ssp245[thres].groupby(df_ssp245.time.dt.year).sum()
    
    #Smooth
    window = np.ones(10) / 10
    
    df_lis_ssp245=pd.DataFrame(columns=["year",thres])
    df_lis_ssp245["year"]=df_year_ssp245.index.values[5:-5]
    df_lis_ssp245[thres]=np.convolve(df_year_ssp245.values,window,mode='same')[5:-5]
    
    
    df_ssp245.to_csv(dirres+"/"+station+"_"+modele+"_binary_over3days_NOAA_"+thres+"_HIX_ssp245_1990-2100.csv")
    df_lis_ssp245.to_csv(dirres+"/"+station+"_"+modele+"_nb_days_over3days_NOAA_"+thres+"_HIX_ssp245_1990-2100.csv")

    
    
    
    times_ssp585=data_ssp585.time.values
    
    i=0
    
    j=0

    
    while i<(times_ssp585.shape[0]):
        
        
        if s_ssp585[i]==1:
            
            j+=1
            
        
        else:
            
            if j<3 and j>0 and i>2:
                
                s_ssp585[i-1]=0
                s_ssp585[i-2]=0
            
            j=0
        
            
        i+=1
    
    df_ssp585=pd.DataFrame(columns=["time",thres])
    
    df_ssp585['time']=pd.to_datetime(data_ssp585.time)
    df_ssp585[thres]=s_ssp585
    
    
    df_year_ssp585=df_ssp585[thres].groupby(df_ssp585.time.dt.year).sum()
    
    #Smooth
    window = np.ones(10) / 10
    
    df_lis_ssp585=pd.DataFrame(columns=["year",thres])
    df_lis_ssp585["year"]=df_year_ssp585.index.values[5:-5]
    df_lis_ssp585[thres]=np.convolve(df_year_ssp585.values,window,mode='same')[5:-5]
    
    df_ssp585.to_csv(dirres+"/"+station+"_"+modele+"_binary_over3days_NOAA_"+thres+"_HIX_ssp585_1990-2100.csv")
    df_lis_ssp585.to_csv(dirres+"/"+station+"_"+modele+"_nb_days_over3days_NOAA_"+thres+"_HIX_ssp585_1990-2100.csv")
    


# In[5]:


for station in Stations:
    
    for modele in Modeles:

        # files with: 1 if HIX>=thres, 0 otherwise
        binary_ssp245=pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_NOAA_thresholds_HIX_ssp245_1990-2100.csv")
        binary_ssp585=pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_NOAA_thresholds_HIX_ssp585_1990-2100.csv")
    
        for thres in ['s1','s2','s3']:
        
            count_periods_over3days(binary_ssp245,binary_ssp585,thres)








####################################################################################
# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import pandas as pd
import matplotlib.ticker as ticker
from tqdm import tqdm


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
    for i in tqdm(range(n)):
        res[i,:]=X[randind[i],:]
    return res


# In[5]:


def calcul_nb_days_heatwaves(station,scenario):
    
    
        
    '''
    This method aims at giving uncertainties and statistics by bootstrap over number of days under HHW1 conditions
    We test different thresholds:
    s1 : CAUTION
    s2 : DANGER
    s3 : EXTREME DANGER
    '''
    
    data_stations_s1=np.empty((len(Modeles),len(range_date)))
    data_stations_s2=np.empty((len(Modeles),len(range_date)))
    data_stations_s3=np.empty((len(Modeles),len(range_date)))

    
    i=0
    
    
    for modele in Modeles:
        
        s1_data=pd.read_csv(dirres+"/"+station+"_"+modele+"_nb_days_over3days_NOAA_s1_HIX_"+scenario+"_1990-2100.csv",header=0,names=['year','s1'])
        s2_data=pd.read_csv(dirres+"/"+station+"_"+modele+"_nb_days_over3days_NOAA_s2_HIX_"+scenario+"_1990-2100.csv",header=0,names=['year','s2'])
        s3_data=pd.read_csv(dirres+"/"+station+"_"+modele+"_nb_days_over3days_NOAA_s3_HIX_"+scenario+"_1990-2100.csv",header=0,names=['year','s3'])

        data_stations_s1[i,:]=s1_data.s1.values[:]
        data_stations_s2[i,:]=s2_data.s2.values[:]
        data_stations_s3[i,:]=s3_data.s3.values[:]

        
        i+=1
        
    n=1000000
    
    
    samples_s1=bootstrap(data_stations_s1,n)
    samples_s2=bootstrap(data_stations_s2,n)
    samples_s3=bootstrap(data_stations_s3,n)

    
    q10 = np.nanquantile(samples_s1,0.10,axis=0)
    q90 = np.nanquantile(samples_s1,0.90,axis=0)
    median= np.nanmedian(samples_s1,axis=0)
    
    q10_df=pd.DataFrame(index=range_date.tolist())
    q10_df['nb_days']=q10
    q90_df=pd.DataFrame(index=range_date.tolist())
    q90_df['nb_days']=q90
    median_df=pd.DataFrame(index=range_date.tolist())
    median_df['nb_days']=median
    
    q10_df.to_csv(dirres+station+"_q10_s1_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.csv")
    q90_df.to_csv(dirres+station+"_q90_s1_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.csv")
    median_df.to_csv(dirres+station+"_q50_s1_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.csv")
    np.save(dirres+station+"_samples_s1_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.npy",samples_s1)
    
    q10 = np.nanquantile(samples_s2,0.10,axis=0)
    q90 = np.nanquantile(samples_s2,0.90,axis=0)
    median= np.nanmedian(samples_s2,axis=0)
    
    q10_df=pd.DataFrame(index=range_date.tolist())
    q10_df['nb_days']=q10
    q90_df=pd.DataFrame(index=range_date.tolist())
    q90_df['nb_days']=q90
    median_df=pd.DataFrame(index=range_date.tolist())
    median_df['nb_days']=median
    
    q10_df.to_csv(dirres+station+"_q10_s2_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.csv")
    q90_df.to_csv(dirres+station+"_q90_s2_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.csv")
    median_df.to_csv(dirres+station+"_q50_s2_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.csv")
    np.save(dirres+station+"_samples_s2_days_NOAA_over3days_thresholds_HIX_"+scenario+"_CMIP6.npy",samples_s2)
    
    q10 = np.nanquantile(samples_s3,0.10,axis=0)
    q90 = np.nanquantile(samples_s3,0.90,axis=0)
    median= np.nanmedian(samples_s3,axis=0)
    
    q10_df=pd.DataFrame(index=range_date.tolist())
    q10_df['nb_days']=q10
    q90_df=pd.DataFrame(index=range_date.tolist())
    q90_df['nb_days']=q90
    median_df=pd.DataFrame(index=range_date.tolist())
    median_df['nb_days']=median
    
    q10_df.to_csv(dirres+station+"_q10_s3_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.csv")
    q90_df.to_csv(dirres+station+"_q90_s3_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.csv")
    median_df.to_csv(dirres+station+"_q50_s3_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.csv")
    np.save(dirres+station+"_samples_s3_days_over3days_NOAA_thresholds_HIX_"+scenario+"_CMIP6.npy",samples_s3)
    


# In[6]:


for station in Stations:
    for scenario in Scenarios:
        calcul_nb_days_heatwaves(station,scenario)


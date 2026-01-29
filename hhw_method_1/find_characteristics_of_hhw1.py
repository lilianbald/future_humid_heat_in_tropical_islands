#!/usr/bin/env python
# coding: utf-8


# In[1]:


import numpy as np
import pandas as pd
from datetime import datetime 


# In[2]:


dirres= "/CMIP6/heatwaves/"
dirout = "/CMIP6/heatwaves/"


# In[3]:


Stations=[
    "Gillot","Pamandzi","Tromelin"
    ]



# In[1]:


# Here we deduce max intensity and duration of each HHW1


# In[4]:


for station in Stations:
    
    Modeles=["ACCESS-CM2","CanESM5",
         "CNRM-ESM2-1",
        "CNRM-CM6-1",
        "CNRM-CM6-1-HR",
         "EC-Earth3","EC-Earth3-CC","EC-Earth3-Veg-LR",
         "IPSL-CM6A-LR","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR",
         "MPI-ESM1-2-LR",
         "NorESM2-MM"
        ]
    
    for modele in Modeles:
        
        data_ssp245 = pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_over3days_NOAA_s2_HIX_ssp245_1990-2100.csv")
        data_ssp585= pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_over3days_NOAA_s2_HIX_ssp585_1990-2100.csv")
        HIX_ssp245=pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_NOAA_thresholds_HIX_ssp245_1990-2100.csv")
        HIX_ssp585=pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_NOAA_thresholds_HIX_ssp585_1990-2100.csv")
        
        
        data_ssp245["time"]=pd.to_datetime(data_ssp245.time,errors='coerce')
        data_ssp585["time"]=pd.to_datetime(data_ssp585.time,errors='coerce')
        
        HIX_ssp245["time"]=pd.to_datetime(HIX_ssp245.time,errors='coerce')
        HIX_ssp585["time"]=pd.to_datetime(HIX_ssp585.time,errors='coerce')
        
        year_ssp245=[]
        lenght_ssp245=[]
        intensity_ssp245=[]
        
        day=data_ssp245["time"].values[0]
        
        if modele == "CNRM-CM6-1-HR":
            
            date_fin="2100-12-31T10:30:00"
            
        else:
            
            date_fin="2100-12-31T12:00:00"
        
        while not day == pd.to_datetime(date_fin):
            
            #print(day)
            
            if data_ssp245.loc[data_ssp245["time"]==day].s2.values==1:
                
                year_ssp245.append(pd.to_datetime(day,errors="coerce").year)
                
                lenght_hhw=0
                
                intensity_hhw_list=[]
                                
                while (data_ssp245.loc[data_ssp245["time"]==day].s2.values==1) and not (day == pd.to_datetime(date_fin)) :
                    
                    lenght_hhw+=1
                                        
                    intensity_hhw_list.append(HIX_ssp245.loc[HIX_ssp245["time"]==day].HI.values)
                    
                    iday=data_ssp245.loc[data_ssp245["time"]==day].index
                    
                    
                    if not iday>=data_ssp245.time.values.shape[0]-1:
                    
                        day=data_ssp245.iloc[iday+1].time.values[0]
                    
                    
                lenght_ssp245.append(lenght_hhw)
                intensity_ssp245.append(np.max(intensity_hhw_list))
                        
            iday=data_ssp245.loc[data_ssp245["time"]==day].index
            
            if not iday>=data_ssp245.time.values.shape[0]-1:
                
                day=data_ssp245.iloc[iday+1].time.values[0]
        
        if lenght_ssp245[-1]<3:
            
            year_ssp245=year_ssp245[:-1]
            lenght_ssp245=lenght_ssp245[:-1]
            intensity_ssp245=intensity_ssp245[:-1]

        
        df_ssp245=pd.DataFrame()
        
        df_ssp245["year"]=np.array(year_ssp245).astype(int)
        df_ssp245["length"]=np.array(lenght_ssp245).astype(int)
        df_ssp245["intensity"]=intensity_ssp245
        
        print(station)
        print(modele)
        print('ssp245 done')
        df_ssp245.to_csv(dirout+station+"_"+modele+"_carac_hhw_NOAA_s2_over3days_ssp245_1990-2100.csv")
        
        
        
        
        
        year_ssp585=[]
        lenght_ssp585=[]
        intensity_ssp585=[]
        
        day=data_ssp585["time"].values[0]
        
        while not day == pd.to_datetime(date_fin):
            
            if data_ssp585.loc[data_ssp585["time"]==day].s2.values==1:
                
                year_ssp585.append(pd.to_datetime(day,errors="coerce").year)
                
                lenght_hhw=0
                
                intensity_hhw_list=[]
                                
                while (data_ssp585.loc[data_ssp585["time"]==day].s2.values==1) and not (day == pd.to_datetime(date_fin)) :
                    
                    lenght_hhw+=1
                                        
                    intensity_hhw_list.append(HIX_ssp585.loc[HIX_ssp585["time"]==day].HI.values)
                    
                    iday=data_ssp585.loc[data_ssp585["time"]==day].index
                    
                    
                    if not iday>=data_ssp585.time.values.shape[0]-1:
                    
                        day=data_ssp585.iloc[iday+1].time.values[0]
                    
                    
                lenght_ssp585.append(lenght_hhw)
                intensity_ssp585.append(np.max(intensity_hhw_list))
                        
            iday=data_ssp585.loc[data_ssp585["time"]==day].index
            
            if not iday>=data_ssp585.time.values.shape[0]-1:
                
                day=data_ssp585.iloc[iday+1].time.values[0]
        
        if lenght_ssp585[-1]<3:
            
            year_ssp585=year_ssp585[:-1]
            lenght_ssp585=lenght_ssp585[:-1]
            intensity_ssp585=intensity_ssp585[:-1]


        df_ssp585=pd.DataFrame()
        
        df_ssp585["year"]=np.array(year_ssp585).astype(int)
        df_ssp585["length"]=np.array(lenght_ssp585).astype(int)
        df_ssp585["intensity"]=intensity_ssp585
        
        print(station)
        print(modele)
        print('ssp585 done')
        
        df_ssp585.to_csv(dirout+station+"_"+modele+"_carac_hhw_NOAA_s2_over3days_ssp585_1990-2100.csv")
        


# In[2]:


# For calendar reason, KACE-1-0-G is processed hereafter


# In[5]:


for station in Stations:
    
    Modeles=[
         "KACE-1-0-G"
        ]
    
    for modele in Modeles:
        
        data_ssp245 = pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_over3days_NOAA_s2_HIX_ssp245_1990-2100.csv")
        data_ssp585= pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_over3days_NOAA_s2_HIX_ssp585_1990-2100.csv")
        HIX_ssp245=pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_NOAA_thresholds_HIX_ssp245_1990-2100.csv")
        HIX_ssp585=pd.read_csv(dirres+"/"+station+"_"+modele+"_binary_NOAA_thresholds_HIX_ssp585_1990-2100.csv")
        
        data_ssp245["time"]=pd.to_datetime(data_ssp245.time,errors='coerce')
        data_ssp585["time"]=pd.to_datetime(data_ssp585.time,errors='coerce')
        
        HIX_ssp245["time"]=pd.to_datetime(HIX_ssp245.time,errors='coerce')
        HIX_ssp585["time"]=pd.to_datetime(HIX_ssp585.time,errors='coerce')
        
        
        #remove days 31
        index_31=data_ssp245.loc[data_ssp245.time.dt.day==31].index.values
        
        data_ssp245.drop(index=index_31,inplace=True)
        data_ssp245.set_index(np.arange(data_ssp245.time.values.shape[0]),inplace=True)        
        data_ssp585.drop(index=index_31,inplace=True)
        data_ssp585.set_index(np.arange(data_ssp585.time.values.shape[0]),inplace=True)
        HIX_ssp245.drop(index=index_31,inplace=True)
        HIX_ssp245.set_index(np.arange(HIX_ssp245.time.values.shape[0]),inplace=True)
        HIX_ssp585.drop(index=index_31,inplace=True)
        HIX_ssp585.set_index(np.arange(HIX_ssp585.time.values.shape[0]),inplace=True)

        
        year_ssp245=[]
        lenght_ssp245=[]
        intensity_ssp245=[]
        
        day=data_ssp245["time"].values[0]
        
        while not day == pd.to_datetime("2100-12-30T00:00:00"):
                        
            if data_ssp245.loc[data_ssp245["time"]==day].s2.values==1:
                
                year_ssp245.append(pd.to_datetime(day,errors="coerce").year)
                
                lenght_hhw=0
                
                intensity_hhw_list=[]
                                
                while (data_ssp245.loc[data_ssp245["time"]==day].s2.values==1) and not (day == pd.to_datetime("2100-12-30T00:00:00")) :
                    
                    lenght_hhw+=1
                                        
                    intensity_hhw_list.append(HIX_ssp245.loc[HIX_ssp245["time"]==day].HI.values)
                    
                    iday=data_ssp245.loc[data_ssp245["time"]==day].index
                    
                    
                    if not iday>=data_ssp245.time.values.shape[0]-1:
                    
                        day=data_ssp245.iloc[iday+1].time.values[0]
                    
                    
                lenght_ssp245.append(lenght_hhw)
                intensity_ssp245.append(np.max(intensity_hhw_list))
            
            #print(day)
            iday=data_ssp245.loc[data_ssp245["time"]==day].index
                        
            if not iday>=data_ssp245.time.values.shape[0]-1:
                
                day=data_ssp245.iloc[iday+1].time.values[0]
        
        if lenght_ssp245[-1]<3:# on ne compte pas les vagues de chaleur de moins de 3 jours --> pb à l'approche de 31/12/2100
            
            year_ssp245=year_ssp245[:-1]
            lenght_ssp245=lenght_ssp245[:-1]
            intensity_ssp245=intensity_ssp245[:-1]

        
        df_ssp245=pd.DataFrame()
        
        df_ssp245["year"]=np.array(year_ssp245).astype(int)
        df_ssp245["length"]=np.array(lenght_ssp245).astype(int)
        df_ssp245["intensity"]=intensity_ssp245
        
        #print(df_ssp245)
        
        print(station)
        print(modele)
        print('ssp245 done')
        df_ssp245.to_csv(dirout+station+"_"+modele+"_carac_hhw_NOAA_s2_over3days_ssp245_1990-2100.csv")
        
        
        
        
        
        year_ssp585=[]
        lenght_ssp585=[]
        intensity_ssp585=[]
        
        day=data_ssp585["time"].values[0]
        
        while not day == pd.to_datetime("2100-12-30T00:00:00"):
            
            if data_ssp585.loc[data_ssp585["time"]==day].s2.values==1:
                
                year_ssp585.append(pd.to_datetime(day,errors="coerce").year)
                
                lenght_hhw=0
                
                intensity_hhw_list=[]
                                
                while (data_ssp585.loc[data_ssp585["time"]==day].s2.values==1) and not (day == pd.to_datetime("2100-12-30T00:00:00")) :
                    
                    lenght_hhw+=1
                                        
                    intensity_hhw_list.append(HIX_ssp585.loc[HIX_ssp585["time"]==day].HI.values)
                    
                    iday=data_ssp585.loc[data_ssp585["time"]==day].index
                    
                    
                    if not iday>=data_ssp585.time.values.shape[0]-1:
                    
                        day=data_ssp585.iloc[iday+1].time.values[0]
                    
                    
                lenght_ssp585.append(lenght_hhw)
                intensity_ssp585.append(np.max(intensity_hhw_list))
                        
            iday=data_ssp585.loc[data_ssp585["time"]==day].index
            
            if not iday>=data_ssp585.time.values.shape[0]-1:
                
                day=data_ssp585.iloc[iday+1].time.values[0]
        
        if lenght_ssp585[-1]<3:# on ne compte pas les vagues de chaleur de moins de 3 jours --> pb à l'approche de 31/12/2100
            
            year_ssp585=year_ssp585[:-1]
            lenght_ssp585=lenght_ssp585[:-1]
            intensity_ssp585=intensity_ssp585[:-1]


        df_ssp585=pd.DataFrame()
        
        df_ssp585["year"]=np.array(year_ssp585).astype(int)
        df_ssp585["length"]=np.array(lenght_ssp585).astype(int)
        df_ssp585["intensity"]=intensity_ssp585
        
        print(station)
        print(modele)
        print('ssp585 done')        
        
        df_ssp585.to_csv(dirout+station+"_"+modele+"_carac_hhw_NOAA_s2_over3days_ssp585_1990-2100.csv")
        

###########################################################################################

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


dirin       = '/CMIP6/heatwaves/'
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


Stations=["Gillot",
    "Pamandzi","Tromelin"
    ]


# In[ ]:


# Here we juste aggregate CMIP6 data of HHW1s, sorted by year. We keep the year, the intensity and the duration of the event.


# In[4]:


for station in Stations:
    
    print(station+'\n')
    
    df_station_ssp245=pd.DataFrame(columns=["year","length","intensity"])
    df_station_ssp585=pd.DataFrame(columns=["year","length","intensity"])

    
    for modele in Modeles:
        
        print(modele)
        
        data_ssp245 = pd.read_csv(dirres+station+'_'+modele+'_carac_hhw_NOAA_s2_over3days_ssp245_1990-2100.csv')
        data_ssp585 = pd.read_csv(dirres+station+'_'+modele+'_carac_hhw_NOAA_s2_over3days_ssp585_1990-2100.csv')
        
        print(data_ssp245)
                
        df_ssp245=pd.DataFrame()
        df_ssp585=pd.DataFrame()
        
        if not data_ssp245.year.values.shape[0]==0:
            
            df_ssp245["year"]=data_ssp245.year.values.astype(int)
            df_ssp245["length"]=data_ssp245.length.values.astype(int)
            df_ssp245["intensity"]=data_ssp245.intensity.values
            
            df_ssp245.to_csv(dirres+"/"+station+"_"+modele+"_carac_heatwaves_DATAFRAME_HIX_ssp245_1990-2100.csv")
            
            df_station_ssp245=pd.concat([df_station_ssp245,df_ssp245])


        if not data_ssp585.year.values.shape[0]==0:
            
            df_ssp585["year"]=data_ssp585.year.values.astype(int)
            df_ssp585["length"]=data_ssp585.length.values.astype(int)
            df_ssp585["intensity"]=data_ssp585.intensity.values
                
            df_ssp585.to_csv(dirres+"/"+station+"_"+modele+"_carac_heatwaves_DATAFRAME_HIX_ssp585_1990-2100.csv")
        
            df_station_ssp585=pd.concat([df_station_ssp585,df_ssp585])
        
            
    df_station_ssp245.to_csv(dirres+station+"_ALL_MODELES_carac_hhw_NOAA_s2_over3days_HIX_ssp245_1990-2100.csv")
    df_station_ssp585.to_csv(dirres+station+"_ALL_MODELES_carac_hhw_NOAA_s2_over3days_HIX_ssp585_1990-2100.csv") 


############################################################################################

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


def calcul_carac_heatwaves(station,scenario):
    
    '''
    This method aims at getting characteristics of HHW1 (max intensity and duration) from all models together
    and for each decadal in period 1990-2100.
    '''

    
    i=0
    
        
    data_station=pd.read_csv(dirres+"/"+station+"_ALL_MODELES_carac_hhw_NOAA_s2_over3days_HIX_"+scenario+"_1990-2100.csv",header=0,names=['year','length','intensity'])
    
    #print(data_station)
    
    data_station_sort=data_station.sort_values(by='year')
    
    dec=[]
    
    for year in data_station_sort.year.values:
        
        if not year == 1990:
            
            dec_year=((year-1)//10)*10 + 5
        
            dec.append(dec_year)
            
        else:
            
            dec.append(1995)
            
    
    data_station_sort["dec"]=dec
    
    
    for dec_value in np.unique(data_station_sort.dec.values):
        
        
        data_dec=data_station_sort[data_station_sort.dec==dec_value]
        
        np.save(dirres+station+"_length_heatwaves_NOAA_s2_over3days_HIX_"+scenario+"_CMIP6_"+str(dec_value)+".npy",data_dec.length.values)   
        np.save(dirres+station+"_intensity_heatwaves_NOAA_s2_over3days_HIX_"+scenario+"_CMIP6_"+str(dec_value)+".npy",data_dec.intensity.values)    


# In[5]:


for station in Stations:
    for scenario in Scenarios:
        calcul_carac_heatwaves(station,scenario)


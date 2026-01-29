#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import netCDF4 as nc
import cftime
import pandas as pd
from datetime import timedelta
from datetime import datetime 


# In[2]:


period1     = "1990-2014"
period2     = "2015-2100" 
dirin       = '/OBS/HIM/'
dirres= "/OBS/heatwaves/"


date_format = "%Y-%m-%d"


# In[3]:


# Definition of the function to find heat waves
# This version uses thresholds S1, S2, and S3 to detect the start, continuation, and end of heat waves. 
# Days when the temperature is above S2 are added to the heatwave, and the intensity of the heatwave is calculated
#as the sum of the maximum HIs above S2. 
#The end of a heatwave is detected when the HI falls below


# In[4]:


def find_peaks(Spic,tas, dates):
    '''
    Find the dates of heat peaks in a series of average HI data.
    Spic is the threshold above which a temperature is considered high enough to detect the start of a heatwave.
    '''

    heatwaves_dates=[dates[i] for i in range(len(tas)) if tas[i]>=Spic]
    return filter_first_dates(heatwaves_dates)   # eliminer les doublons de pics de chaleur depassant S1 à des dates successives


# In[5]:


# Definition of the function to eliminate duplicate heat peaks
def filter_first_dates(dates):
    """
    in a table of daily dates, keep only the first date of each consecutive sequence, separated by a single day    """
    filtered_dates = [dates[0]] 
    for i in range(1, len(dates)):
        if (pd.to_datetime(dates[i]) - pd.to_datetime(dates[i-1])).days > 1:  
   
            filtered_dates.append(dates[i]) 


    return filtered_dates


# In[6]:


def find_heatwave_deb(Sint,Sdeb,Spic,dates_pic,tas,dates):
    """
    determines the start date of a heatwave in a series of average HI data (HI).
    Spic is the threshold above which a HI is considered high enough to detect the start of a heatwave.
    Sdeb is the threshold above which a HI is considered high enough to add a day to the heatwave.
    Sint is the threshold below which a HI is considered relatively cold for detecting the end of a heatwave.
    """

    dates_deb=[]
    list_dates = dates.tolist()
    

    for i in range(len(dates_pic)):
        booleen = True
    
        
        d       = dates_pic[i]
        J       = list_dates.index(d)
        
        print("date du pic ", d)
        print("J ",J)

        while booleen == True:
            Jmoins1 = J-1
           

            if tas[Jmoins1] < Sint:
                datedeb = dates.values[J]
                booleen=False
            elif (tas[Jmoins1] >= Sint and tas[Jmoins1] < Sdeb):
                Jmoins2 = J - 2
                Jmoins3 = J - 3
        
                if (tas[Jmoins2] < Sdeb and tas[Jmoins3] < Sdeb) or tas[Jmoins2] < Sint:
                    datedeb = dates.values[J]
                    booleen=False
                else:
                    J = J - 1  
            else:   
              
                J = J -1

        dates_deb.append((datedeb))

    return filter_first_dates(dates_deb)


# In[7]:


def find_heatwave_fin(Sint,Sdeb,Spic,dates_pic,tas,dates):
    """
    determines the start date of a heatwave in a series of average HI data (HI).
    Spic is the threshold above which a HI is considered high enough to detect the start of a heatwave.
    Sdeb is the threshold above which a HI is considered high enough to add a day to the heatwave.
    Sint is the threshold below which a HI is considered relatively cold for detecting the end of a heatwave.
    """

    dates_fin=[]
    list_dates = dates.tolist() 

    for i in range(len(dates_pic)):
        booleen = True
        d       = dates_pic[i]
        J       = list_dates.index(d)
     
        if J+3 < len(dates):
            while ((booleen == True) and (J+3 < len(dates))):
                Jplus1 = J+1
          
                if Jplus1 > len(dates) - 2 : 
                    tas[Jplus1] = 0

                if tas[Jplus1] < Sint:
                    datefin = dates.values[J]
                    booleen=False
                elif (tas[Jplus1] >= Sint and tas[Jplus1] < Sdeb):
                    Jplus2 = J + 2
                    Jplus3 = J + 3
                    if Jplus2 > len(dates) - 2 : 
                        tas[Jplus2] = 0
                    if Jplus3 > len(dates) - 2 :
                        tas[Jplus3] = 0
        
                    if (tas[Jplus2] < Sdeb and tas[Jplus3] < Sdeb) or tas[Jplus2] < Sint:
                        datefin = dates.values[J]
                        booleen=False
                    else: 
                        J = J + 1  
                else:   
                    J = J + 1

            dates_fin.append((datefin))

    return filter_first_dates(dates_fin)  


# In[8]:


def heatwaves(Sdeb, Spic, dates_debut,dates_fin,tas,dates):
    """
    Characterises the duration and intensity of humid heat waves in a series of maximum HI data.
    """

    # initialisations
    num_rows = 0 
    parameters = []
    years_tab = [] 
    for date in dates:
        annee = date.year
        if annee not in years_tab:
            years_tab.append(annee)
    cumul_annuel = [0] * len(years_tab)  #  List for storing annual heatwave day totals
    cumul_annuel_hw=[0] * len(years_tab)

    if (dates_debut[-1]>dates_fin[-1]):
        dates_fin.append(datetime(2014, 12, 31))
        
                
    if len(dates_debut)>len(dates_fin):
        dates_debut=dates_debut[:-1]
 
    for i in range(len(dates_debut)):
        debut = pd.to_datetime(dates_debut[i])
        Jdeb = np.where(dates == debut)[0] 
       
       
        fin = pd.to_datetime(dates_fin[i])
        Jfin = np.where(dates == fin)[0] 
        
   
        if debut.year != fin.year:
            if fin.year - debut.year >1:
                
                fin_annee_civile = datetime(debut.year, 12, 31)
                duree1 = (fin_annee_civile - debut).days + 1  
                cumul_annuel[debut.year - years_tab[0]] += duree1

                if duree>=3:
                    
                    cumul_annuel_hw[debut.year - years_tab[0]] += duree1

                debut_annee_suivante = datetime(fin.year, 1, 1)
                    
                duree2 = (fin - debut_annee_suivante).days + 1 
                
                cumul_annuel[fin.year - years_tab[0]] += duree2
                if duree2>=3:
                    cumul_annuel_hw[fin.year - years_tab[0]] += duree2
                    
                hw_year=debut.year+1
                
                while not (hw_year == fin.year):
                    cumul_annuel[hw_year- years_tab[0]]+=365
                    cumul_annuel_hw[hw_year- years_tab[0]]+=365
            
                    hw_year+=1
                    
            else:
                fin_annee_civile = datetime(debut.year, 12, 31)
                    
                duree1 = (fin_annee_civile - debut).days + 1  
                cumul_annuel[debut.year - years_tab[0]] += duree1
                
                if duree1>=3:
                    cumul_annuel_hw[debut.year - years_tab[0]] += duree1

              
                debut_annee_suivante = datetime(fin.year, 1, 1)
                    
                duree2 = (fin - debut_annee_suivante).days + 1  
                cumul_annuel[fin.year - years_tab[0]] += duree2
                
                if duree2>=3:
                    cumul_annuel_hw[fin.year - years_tab[0]] += duree2
        else:
            duree = (fin - debut).days + 1 
            cumul_annuel[debut.year - years_tab[0]] += duree
            
            if duree>=3:
                cumul_annuel_hw[debut.year - years_tab[0]] += duree

        duree = Jfin[0] - Jdeb[0] + 1 
        if duree >= 3: 
            intensite_max = round( np.amax(tas[ Jdeb[0]:Jfin[0]+1 ]), 2 )
       
            ecarts = sum((tas - Sdeb) for tas in tas[ Jdeb[0]:Jfin[0]+1 ] if tas >= Sdeb )  
            severite = round( ecarts / (Spic - Sdeb ), 2)
       
            annee_episode = pd.to_datetime(dates_debut[i]).year
            parameters.append((annee_episode,duree,intensite_max,severite))
            num_rows += 1
    
    return parameters, num_rows, cumul_annuel, cumul_annuel_hw


# In[9]:


def calcul_heatwaves(station):
    
    ficin1_path   = os.path.join(dirin,station+"_obs_1985-2014_HIM.csv")
    print('History input file is reading from ' + ficin1_path)  
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    hist_nc = pd.read_csv(ficin1_path)

    T_hist = hist_nc.loc[pd.to_datetime(hist_nc['time'],format="%Y-%m-%d") >= "1990-01-01"][varname][:]
 
    
    time_hist_vals = hist_nc.loc[pd.to_datetime(hist_nc['time'],format="%Y-%m-%d") >= "1990-01-01"][vartime]

    
    hist_dates = pd.to_datetime(time_hist_vals)
  
    annee_debut_hist= hist_dates.dt.year.values[0]  #First year
    annee_fin_hist = hist_dates.dt.year.values[-1]  #Last year


    # --- Define thresholds from historical data
    #S1 is the threshold above which a HI is considered high enough to detect the start of a heatwave.
    #S2 is the threshold above which a HI is considered high enough to add a day to the heatwave.
    #S3 is the threshold below which a HI is considered relatively cool enough to detect the end of a heatwave.
    s1 = np.nanpercentile(T_hist.values, 99.5)
    s2 = np.nanpercentile(T_hist.values, 97.5)
    s3 = np.nanpercentile(T_hist.values, 95.0)

    # Historical observed HHW2
 
    
    pics_hist_all = find_peaks(s1,T_hist.values,hist_dates.values)

    pics_hist     = filter_first_dates(pics_hist_all)

    formatted_pics_hist = pics_hist

    print("-----------HI of peaks in HIST")
    mask_pics = np.in1d(hist_dates.values, pics_hist)
    pics_values = T_hist[mask_pics]

    

    #print("----------DATES Beginning EPISODES HIST")
    dates_deb_hist=find_heatwave_deb(s3,s2,s1,pics_hist,T_hist.values,hist_dates)
    formatted_dates_deb_hist = dates_deb_hist
    #print(formatted_dates_deb_hist)

    #print("----------DATES end EPISODES HIST")
    dates_fin_hist=find_heatwave_fin(s3,s2,s1,pics_hist,T_hist.values,hist_dates)
    formatted_dates_fin_hist = dates_fin_hist

    #print("----------characteristics HHW2 hist obs")
    heatwaves_hist, num_vagues_hist, cumul_annuel_hist, cumul_annuel_hist_hw = heatwaves(s2,s1,dates_deb_hist,dates_fin_hist,T_hist,hist_dates)
    #print(heatwaves_hist)
    print("nbr of HHW2", num_vagues_hist)
    print("days/year of HHW2", cumul_annuel_hist_hw)

    
    #Saving
    np.save(dirres+"/"+station+"_heatwaves_HIM_obs.npy", heatwaves_hist)
    np.save(dirres+"/"+station+"_num_vagues_HIM_obs.npy", num_vagues_hist)
    np.save(dirres+"/"+station+"_annee_debut_HIM_obs.npy", annee_debut_hist)
    np.save(dirres+"/"+station+"_annee_fin_HIM_obs.npy", annee_fin_hist)
    np.save(dirres+"/"+station+"_cumul_annuel_HIM_obs.npy", cumul_annuel_hist)
    np.save(dirres+"/"+station+"_cumul_annuel_hw_HIM_obs.npy", cumul_annuel_hist_hw)

    #return thresholds to be stored
    return s1,s2,s3


# In[10]:


Stations=["Gillot",
    "Pamandzi","Tromelin"
    ]

S1=pd.DataFrame()
S2=pd.DataFrame()
S3=pd.DataFrame()

for station in Stations:
    
    q995=[]
    q975=[]
    q950=[]
    
    
    s1,s2,s3=calcul_heatwaves(station)
        
    q995.append(s1)
    q975.append(s2)
    q950.append(s3)
        
    S1[station]=q995
    S2[station]=q975
    S3[station]=q950

#saving thresolds
S1.to_csv(dirres+"/seuils/q995.csv")
S2.to_csv(dirres+"/seuils/q975.csv")
S3.to_csv(dirres+"/seuils/q950.csv")


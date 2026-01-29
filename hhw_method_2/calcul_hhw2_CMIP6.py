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
dirin       = '/CMIP6/HIM/'
dirres= "/CMIP6/heatwaves/"
vartime= 'time'
list_scenarios = ["ssp245","ssp585"]
varname     = "HI"

date_format = "%Y-%m-%d"


# In[3]:


# Definition of the function to find heat waves
# This version uses thresholds S1, S2, and S3 to detect the start, continuation, and end of heat waves. 
# Days when the temperature is above S2 are added to the heatwave, and the intensity of the heatwave is calculated
#as the sum of the maximum HIs above S2. 
#The end of a heatwave is detected when the HI falls below


# In[4]:


def find_peaks(Spic,tas, dates):
    """
    Find the dates of heat peaks in a series of average HI data.
    Spic is the threshold above which a temperature is considered high enough to detect the start of a heatwave.
    """

    heatwaves_dates=[dates[i] for i in range(len(tas)) if tas[i]>=Spic]
    return filter_first_dates(heatwaves_dates) 


# In[5]:


# Definition de la fonction pour éliminer les doublons de pics de chaleur
def filter_first_dates(dates):
    """
    in a table of daily dates, keep only the first date of each consecutive sequence, separated by a single day   
    """
    filtered_dates = [dates[0]] 
    for i in range(1, len(dates)):
        if (dates[i] - dates[i-1]).days > 1:  
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
         
        if J-3 >= 0 :
            while booleen == True:
                Jmoins1 = J-1
         

                if tas[Jmoins1] < Sint:
                    datedeb = dates[J]
                    booleen=False
                elif (tas[Jmoins1] >= Sint and tas[Jmoins1] < Sdeb):
                    Jmoins2 = J - 2
                    Jmoins3 = J - 3
        
                    if (tas[Jmoins2] < Sdeb and tas[Jmoins3] < Sdeb) or tas[Jmoins2] < Sint:
                        datedeb = dates[J]
                        booleen=False
                    else:
                        J = J - 1  
                else:  
                    J = J -1
            
            dates_deb.append((datedeb))
            
        else:
            
            datedeb=dates[0]

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
         
                    datefin = dates[J]
          
                    booleen=False
                elif (tas[Jplus1] >= Sint and tas[Jplus1] < Sdeb):
         
                    Jplus2 = J + 2
                    Jplus3 = J + 3
                    if Jplus2 > len(dates) - 2 :
                        tas[Jplus2] = 0
                    if Jplus3 > len(dates) - 2 :
                        tas[Jplus3] = 0
        
                    if (tas[Jplus2] < Sdeb and tas[Jplus3] < Sdeb) or tas[Jplus2] < Sint:
                        datefin = dates[J]
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
    cumul_annuel = [0] * len(years_tab)  # List for storing annual heatwave day totals
    cumul_annuel_hw=[0] * len(years_tab)
    
    if (dates_debut[-1]>dates_fin[-1]):
        if dates_fin[-1].year <= 2014:#historique
            dates_fin.append(datetime(2014, 12, 31, 12, 0))
        else:
            
            if modele=="CNRM-CM6-1-HR":
                dates_fin.append(datetime(2100, 12, 31, 10, 30))
            else:
                dates_fin.append(datetime(2100, 12, 31, 12, 0))
                
    if len(dates_debut)>len(dates_fin):
        dates_debut=dates_debut[:-1]
        


    for i in range(len(dates_debut)):
        if (modele == "CanESM5") or (modele == "NorESM2-MM"):
            debut=cftime.DatetimeNoLeap(dates_debut[i].year,dates_debut[i].month,dates_debut[i].day,12,0)
        else:
            debut = dates_debut[i]
        Jdeb = np.where(dates == debut)[0]
      
        if (modele == "CanESM5") or (modele == "NorESM2-MM"):
            fin=cftime.DatetimeNoLeap(dates_fin[i].year,dates_fin[i].month,dates_fin[i].day,12,0)
        else:
            fin = dates_fin[i]
        Jfin = np.where(dates == fin)[0]  
        if debut.year != fin.year:
            if fin.year - debut.year >1:
                
                if (modele == "CanESM5") or (modele == "NorESM2-MM"):
                    fin_annee_civile = cftime.DatetimeNoLeap(debut.year, 12, 31)
                else:
                    fin_annee_civile = datetime(debut.year, 12, 31)
                duree1 = (fin_annee_civile - debut).days + 1 
                cumul_annuel[debut.year - years_tab[0]] += duree1
                
                if duree1>=3:
                    cumul_annuel_hw[debut.year - years_tab[0]] += duree1

               
                if (modele == "CanESM5") or (modele == "NorESM2-MM"):
                    debut_annee_suivante=cftime.DatetimeNoLeap(fin.year, 1, 1)
                else:
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
                if (modele == "CanESM5") or (modele == "NorESM2-MM"):
                    fin_annee_civile = cftime.DatetimeNoLeap(debut.year, 12, 31)
                else:
                    fin_annee_civile = datetime(debut.year, 12, 31)
                    
                duree1 = (fin_annee_civile - debut).days + 1  
                cumul_annuel[debut.year - years_tab[0]] += duree1
                
                if duree1>=3:
                    cumul_annuel_hw[debut.year - years_tab[0]] += duree1

                if (modele == "CanESM5") or (modele == "NorESM2-MM"):
                    debut_annee_suivante=cftime.DatetimeNoLeap(fin.year, 1, 1)
                else:
                    debut_annee_suivante = datetime(fin.year, 1, 1)
                    
                duree2 = (fin - debut_annee_suivante).days + 1 
                cumul_annuel[fin.year - years_tab[0]] += duree2
                
                if duree2>=3:
                    cumul_annuel_hw[fin.year - years_tab[0]] += duree2
        else:
            duree = (fin - debut).days + 1  
            cumul_annuel[debut.year - years_tab[0]] += duree
            
            if duree >=3:
                cumul_annuel_hw[debut.year - years_tab[0]] += duree

        duree = Jfin[0] - Jdeb[0] + 1 
        if duree >= 3:  
            intensite_max = round( np.amax(tas[ Jdeb[0]:Jfin[0]+1 ]), 2 )
            ecarts = sum((tas - Sdeb) for tas in tas[ Jdeb[0]:Jfin[0]+1 ] if tas >= Sdeb )  
            severite = round( ecarts / (Spic - Sdeb ), 2)
            annee_episode = dates_debut[i].year
            parameters.append((annee_episode,duree,intensite_max,severite))
            num_rows += 1
    
    return parameters, num_rows, cumul_annuel,cumul_annuel_hw


# In[9]:


def calcul_heatwaves(modele,station):
    ficin1_path   = os.path.join(dirin,station+"_"+modele+"_BC_3M_HIM_day_1990-2014.nc")
    print('History input file is reading from ' + ficin1_path)  
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    hist_nc = nc.Dataset(ficin1_path)
 
    T_hist = hist_nc.variables[varname][:]
    time_hist_var = hist_nc.variables[vartime]
    
    if modele == "KACE-1-0-G":#different calendar
        time_hist_vals = time_hist_var[:]+0.5
        T_hist=T_hist.filled(np.nan)
    else:
        time_hist_vals = time_hist_var[:]

    time_hist_units = time_hist_var.units
    time_hist_calendar = time_hist_var.calendar
    
    hist_dates = nc.num2date(time_hist_vals, units=time_hist_units, calendar=time_hist_calendar)
        

    hist_nc.close()
    
    annee_debut_hist= hist_dates[0].year 
    annee_fin_hist = hist_dates[-1].year 


   # --- Define thresholds from historical data
    #S1 is the threshold above which a HI is considered high enough to detect the start of a heatwave.
    #S2 is the threshold above which a HI is considered high enough to add a day to the heatwave.
    #S3 is the threshold below which a HI is considered relatively cool enough to detect the end of a heatwave.
    s1 = np.nanpercentile(T_hist, 99.5)
    s2 = np.nanpercentile(T_hist, 97.5)
    s3 = np.nanpercentile(T_hist, 95.0)

    # Historical HHW2
    print("-----------DATES des pics HIST")
    
    
    pics_hist_all = find_peaks(s1,T_hist,hist_dates)

    pics_hist     = filter_first_dates(pics_hist_all)
    
    formatted_pics_hist = [date.strftime(date_format) for date in pics_hist]
  

    mask_pics = np.in1d(hist_dates, pics_hist)
    pics_values = T_hist[mask_pics]



    print("----------DATES beginningEPISODES HIST")
    dates_deb_hist=find_heatwave_deb(s3,s2,s1,pics_hist,T_hist,hist_dates)
    formatted_dates_deb_hist = [date.strftime(date_format) for date in dates_deb_hist]
    print(formatted_dates_deb_hist)

    print("----------DATES end EPISODES HIST")
    dates_fin_hist=find_heatwave_fin(s3,s2,s1,pics_hist,T_hist,hist_dates)
    formatted_dates_fin_hist = [date.strftime(date_format) for date in dates_fin_hist]
    print(formatted_dates_fin_hist)

    print("----------Characteristics HHW2")
    heatwaves_hist, num_vagues_hist, cumul_annuel_hist,cumul_annuel_hist_hw = heatwaves(s2,s1,dates_deb_hist,dates_fin_hist,T_hist,hist_dates)
    print(heatwaves_hist)
 
    print("days/year hist ", cumul_annuel_hist_hw)

    
    #Saving
    np.save(dirres+"/"+station+"_"+modele+"_heatwaves_HIM_hist.npy", heatwaves_hist)
    np.save(dirres+"/"+station+"_"+modele+"_num_vagues_HIM_hist.npy", num_vagues_hist)
    np.save(dirres+"/"+station+"_"+modele+"_annee_debut_HIM_hist.npy", annee_debut_hist)
    np.save(dirres+"/"+station+"_"+modele+"_annee_fin_HIM_hist.npy", annee_fin_hist)
    np.save(dirres+"/"+station+"_"+modele+"_cumul_annuel_HIM_hist.npy", cumul_annuel_hist)
    np.save(dirres+"/"+station+"_"+modele+"_cumul_annuel_hw_HIM_hist.npy", cumul_annuel_hist_hw)


    # FUTURE
    for ind_scenario in range(0,len(list_scenarios)):
        scenario=list_scenarios[ind_scenario]
        print("---------SCENARIO ", scenario, " ----------------")
        ficin2_path   = os.path.join(dirin,station+"_"+modele+"_BC_3M_HIM_day_"+scenario+"_2015-2100.nc")
        print(' input file for the future is reading from ' + ficin2_path)  
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
  
        fut_nc = nc.Dataset(ficin2_path)
        T_fut = fut_nc.variables[varname][:,19,17]
        T_fut = fut_nc.variables[varname][:]
        time_var = fut_nc.variables[vartime]
        time_vals = time_var[:]
        
        if modele == "KACE-1-0-G":#calendar...
            time_vals = time_var[:]+0.5
        else:
            time_vals = time_var[:]
        
        time_units = time_var.units
        time_calendar = time_var.calendar
        
        fut_dates = nc.num2date(time_vals, units=time_units, calendar=time_calendar)

        fut_dates = nc.num2date(time_vals, units=time_units, calendar="proleptic gregorian")
        fut_dates = nc.num2date(time_vals, units=time_units, only_use_cftime_datetimes=True,calendar=time_calendar)
        annee_debut_fut= fut_dates[0].year  # first year
        annee_fin_fut = fut_dates[-1].year  # last year

        years_fut_tab = []
        for date in fut_dates:
            annee = date.year
            if annee not in years_fut_tab:
                years_fut_tab.append(annee)

        fut_nc.close()


        # HHW2 computation
        pics_fut_all = find_peaks(s1,T_fut,fut_dates)
        pics_fut=filter_first_dates(pics_fut_all)
        formatted_pics_fut= [date.strftime(date_format) for date in pics_fut]


        mask_pics = np.in1d(fut_dates, pics_fut)
        pics_values = T_fut[mask_pics]
       
        print("----------DATES Beginning EPISODES FUTURS")
        dates_deb_fut=find_heatwave_deb(s3,s2,s1,pics_fut,T_fut,fut_dates)
        formatted_dates_deb_fut= [date.strftime(date_format) for date in dates_deb_fut ]
        print(formatted_dates_deb_fut)

        print("----------DATES END EPISODES FUTURS")
        dates_fin_fut=find_heatwave_fin(s3,s2,s1,pics_fut,T_fut,fut_dates)
        formatted_dates_fin_fut= [date.strftime(date_format) for date in dates_fin_fut ]
        print(formatted_dates_fin_fut)

        print("----------characteristics HHW2 under scenario ", scenario)
        heatwaves_fut, num_vagues_fut, cumul_annuel_fut,cumul_annuel_fut_hw = heatwaves(s2,s1,dates_deb_fut,dates_fin_fut,T_fut,fut_dates)
        print(heatwaves_fut)
        
        print("nb/year HHW2", cumul_annuel_fut_hw)

        
        
        #Saving
        np.save(dirres+"/"+station+"_"+modele+"_heatwaves_HIM_fut_"+scenario+".npy", heatwaves_fut)
        np.save(dirres+"/"+station+"_"+modele+"_num_vagues_HIM_fut_"+scenario+".npy", num_vagues_fut)
        np.save(dirres+"/"+station+"_"+modele+"_annee_debut_HIM_fut_"+scenario+".npy", annee_debut_fut)
        np.save(dirres+"/"+station+"_"+modele+"_annee_fin_HIM_fut_"+scenario+".npy", annee_fin_fut)
        np.save(dirres+"/"+station+"_"+modele+"_cumul_annuel_HIM_fut_"+scenario+".npy", cumul_annuel_fut)
        np.save(dirres+"/"+station+"_"+modele+"_cumul_annuel_HIM_fut_hw_"+scenario+".npy", cumul_annuel_fut)



    # return thresholds to store it for further analysis    
    return s1,s2,s3


# In[10]:


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
    "Pamandzi","Tromelin"    ]

S1=pd.DataFrame(index=Modeles)
S2=pd.DataFrame(index=Modeles)
S3=pd.DataFrame(index=Modeles)

for station in Stations:
    
    q995=[]
    q975=[]
    q950=[]
    
    for modele in Modeles:
    
        s1,s2,s3=calcul_heatwaves(modele,station)
        
        q995.append(s1)
        q975.append(s2)
        q950.append(s3)
        
    S1[station]=q995
    S2[station]=q975
    S3[station]=q950
    
S1.to_csv(dirres+"/seuils/q995.csv")
S2.to_csv(dirres+"/seuils/q975.csv")
S3.to_csv(dirres+"/seuils/q950.csv")


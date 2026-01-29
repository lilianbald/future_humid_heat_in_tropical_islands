#!/usr/bin/env python
# coding: utf-8

# ## Non-stationnary GEV for return/period level computing

# In[1]:


import pandas as pd
import numpy as np
import time
import xarray as xr
import sys 
import climextremes


# ### Paths

# In[2]:


data_path="./CMIP6/"
covar_path="./CMIP6/covar/"


# ### Modeles, stations, scenarios

# In[4]:


Modeles=["ACCESS-CM2_r2i1p1f1","CanESM5_r1i1p1f1","CNRM-ESM2-1_r1i1p1f2","CNRM-CM6-1_r1i1p1f2","CNRM-CM6-1-HR_r1i1p1f2",
         "EC-Earth3_r1i1p1f1","EC-Earth3-CC_r1i1p1f1","EC-Earth3-Veg-LR_r1i1p1f1","KACE-1-0-G_r1i1p1f1",
         "IPSL-CM6A-LR_r1i1p1f1","MIROC6_r1i1p1f1","MIROC-ES2L_r1i1p1f2","MPI-ESM1-2-HR_r1i1p1f1",
         "MPI-ESM1-2-LR_r1i1p1f1","NorESM2-MM_r1i1p1f1"]

Stations=["Gillot","Pamandzi","Tromelin"]

Scenarios=["ssp245","ssp585"]

loc="Reunion" #South West Indian
#loc='Guyane' for North Atlantic
#loc= ' NC' for Pacific


# ### Data management

# In[3]:


def date_to_datetime(df,time_col):
    '''Converts time data into clean datetime'''
    df[time_col]=pd.to_datetime(df[time_col],errors='coerce')


# In[6]:


def dict_df(var):
    
    '''Create dictionnaries : modele data, covariate
    '''
    
    corr_data=data_path+'/HIX/'
    covar_data=covar_path
    
    scenario_hist='1990-2014'
    
    # create dictionnaries
    
    dict_name='HI_mod'
        
    dict_covar='df_covar'
    
    globals()[dict_name]={}
        
    globals()[dict_covar]={}
    
    for modele_membre in Modeles :
        
        modele=modele_membre[:-9]
        
        globals()[dict_name][modele]={}
                
        globals()[dict_covar][modele]={}
        
        for scenario in Scenarios:
        
        #covariate management
        
            #covar = HI
            globals()[dict_covar][modele][scenario]=pd.read_csv(covar_data+"HI_smoothed_"+scenario+"_"+modele_membre+"_1990-2100_covar_"+loc+".csv",sep=";")
            globals()[dict_covar][modele][scenario].set_index('time',inplace=True)     
            
        
        for station in Stations :
            
            ### corrected data management
            
            globals()[dict_name][modele][station]={}
            
            #future
            
            for scenario in Scenarios:
                            
                globals()[dict_name][modele][station][scenario]=pd.read_csv(corr_data+var+"/"+station+"_"+modele+"_BC_3M_"+var+"_day_"+scenario+"_2015-2100"+".csv",sep=";")
                date_to_datetime(globals()[dict_name][modele][station][scenario],"time")
            
                if (modele == 'CNRM-CM6-1-HR'): # to fix hour issue for this model
                    globals()[dict_name][modele][station][scenario]['time'][1:]=globals()[dict_name][modele][station][scenario]['time'][1:]+pd.Timedelta(1.5,"h")
                
            #hist
            
            globals()[dict_name][modele][station]['hist']=pd.read_csv(corr_data+var+"/"+station+"_"+modele+"_BC_3M_"+var+"_day_"+scenario_hist+".csv",sep=";")
            date_to_datetime(globals()[dict_name][modele][station]['hist'],"time")


# In[7]:


dict_df("HIX") # var = maximum daily HI


# ### Bootstrap to get uncertainties

# In[9]:


def bootstrap(X,Y,n):
    randind=np.random.choice(np.arange(len(X)),(n,len(X)))
    randX=X[randind]
    randY=Y[randind]
    return randind, randX, randY


# ### Methods to compute return levels (RL) and return periods (RP)

# In[15]:


def calcul_extreme_stats_RL(station,period,RP):
    
    '''
    
    Compute return levels of extreme humid heat events
    The gev fit is done with both data from SSP2-4.5 and SSP5-8.5 scenarios,
    considering that gev parameters depend only on covariate values
    
    In the paper: RP = 10 years, period=[1990,2100]
    
    This method returns:
    - Historical averaged return levels
    - Return levels over 1990-2100, in both scenarios SSP2-4.5 and SSP5-8.5
    
    '''
    
    mat_ssp245=np.empty((len(Modeles),(period[1]-period[0])+1))
    mat_ssp245[:]=np.nan
    mat_ssp585=np.empty((len(Modeles),(period[1]-period[0])+1))
    mat_ssp585[:]=np.nan
    
    gev_res={}
    
    mod = 0

    for modele_membre in Modeles :
        
        # modele data
        
        modele=modele_membre[:-9]
        
        data_HI_hist=HI_mod[modele][station]['hist']
        data_HI_fut_ssp245=HI_mod[modele][station]['ssp245']
        data_HI_fut_ssp585=HI_mod[modele][station]['ssp585']
    
        data_HI_hist_fut_ssp245=pd.concat([data_HI_hist,data_HI_fut_ssp245])
        data_HI_hist_fut_ssp585=pd.concat([data_HI_hist,data_HI_fut_ssp585])
        
        # annual maxima

        data_HI_ssp245=data_HI_hist_fut_ssp245.groupby(data_HI_hist_fut_ssp245.time.dt.year).max().HI.loc[period[0]:period[1]]
        data_HI_ssp585=data_HI_hist_fut_ssp585.groupby(data_HI_hist_fut_ssp585.time.dt.year).max().HI.loc[period[0]:period[1]]
    
        data_HI=np.concatenate((data_HI_ssp245.values.reshape(-1,1),data_HI_ssp585.values.reshape(-1,1)),axis=0)
            
        Y=data_HI
    
        # covariate
        
        data_covar_ssp245=df_covar[modele]['ssp245'].HI.loc[period[0]:period[1]]
        data_covar_ssp585=df_covar[modele]['ssp585'].HI.loc[period[0]:period[1]]
        data_covar=np.concatenate((data_covar_ssp245.values.reshape(-1,1),data_covar_ssp585.values.reshape(-1,1)),axis=0)
    
        X=data_covar
    
        
        # GEV fit
        
        gevt=climextremes.fit_gev(Y,X, locationFun = 1, scaleFun = 1, returnPeriod = RP,
                              optimArgs={"method" : "BFGS"},
                              initial={"location" : 40, "scale" : 1, "shape" : -0.01},
                              getParams=True, 
                              #bootSE=True
                             #getFit=True
                             )#likelihood maximisation : algo BFGS (Broyden-Fletcher-Goldfarb-Shanno)
        
        if gevt['info']['failure'][0]==0: #if convergence
            
            return_levels_ssp245=gevt['returnValue'][:(period[1]-period[0])+1]
            return_levels_ssp585=gevt['returnValue'][(period[1]-period[0])+1:]
            
            #storing RL for all models
                
            mat_ssp245[mod,:]=return_levels_ssp245
            mat_ssp585[mod,:]=return_levels_ssp585
                
        mod+=1
    
    # Save results
    
    #RL moyen historique par modèle
    gev_res['RL_hist']={}
    gev_res['RL_hist']['ssp245']=np.nanmean(mat_ssp245[:,:25],axis=1)
    gev_res['RL_hist']['ssp585']=np.nanmean(mat_ssp585[:,:25],axis=1)
        
    #return levels moyen, modèles agrégés
    gev_res['RL_mean']={}
    gev_res['RL_mean']['ssp245']=np.nanmean(mat_ssp245, axis=0)
    gev_res['RL_mean']['ssp585']=np.nanmean(mat_ssp585, axis=0)
    
    data_array_ssp245 = xr.Dataset(dict(
        
        RL =(["year"], gev_res['RL_mean']['ssp245']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RL_mean']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return levels of HI : °C"),

        )
    
    data_array_ssp585 = xr.Dataset(dict(
        
        RL =(["year"], gev_res['RL_mean']['ssp585']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RL_mean']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return levels of HI : °C"),

        )
    
    data_array_ssp245_hist = xr.Dataset(dict(
        
        RL =(["modele"], gev_res['RL_hist']['ssp245']),

        ),
                               
        coords=dict(

        modele=("modele", Modeles),
        
        ),
        
        attrs=dict(units="Return levels of HI : °C"),

        )
    
    data_array_ssp585_hist = xr.Dataset(dict(
        
        RL =(["modele"], gev_res['RL_hist']['ssp585']),

        ),
                               
        coords=dict(

        modele=("modele",Modeles),
        
        ),
        
        attrs=dict(units="Return levels of HI : °C"),

        )
    
    # Saving data
    
    data_array_ssp245_hist.to_netcdf(path=data_path+"GEV/{}y-RL_HI_hist_ssp245_".format(RP)+station+'.nc')
    data_array_ssp585_hist.to_netcdf(path=data_path+"GEV/{}y-RL_HI_hist_ssp585_".format(RP)+station+'.nc')

    data_array_ssp245.to_netcdf(path=data_path+"GEV/{}y-RL_HI_1990-2100_ssp245_".format(RP)+station+'.nc')
    data_array_ssp585.to_netcdf(path=data_path+"GEV/{}y-RL_HI_1990-2100_ssp585_".format(RP)+station+'.nc')


# In[17]:


def calcul_extreme_stats_RP(station,period,RL_hist_ssp245,RL_hist_ssp585,RP):
    
    '''
    
    Compute return periods of historical averaged extreme humid heat events
    The gev fit is done with both data from SSP2-4.5 and SSP5-8.5 scenarios,
    considering that gev parameters depend only on covariate values
    
    In the paper: RP = 10 years, period=[1990,2100], RL_hist_ssp245
    and RL_hist_ssp585 are averaged historical return levels
    
    This method returns:
    - Return periods of historical events over 1990-2100, in both scenarios SSP2-4.5 and SSP5-8.5
    
    '''
    
    mat_ssp245=np.empty((len(Modeles),(period[1]-period[0])+1))
    mat_ssp245[:]=np.nan
    mat_ssp585=np.empty((len(Modeles),(period[1]-period[0])+1))
    mat_ssp585[:]=np.nan
    
    gev_res={}
    
    mod = 0

    for modele_membre in Modeles :
        
        # modele data
                
        modele=modele_membre[:-9]
        
        data_HI_hist=HI_mod[modele][station]['hist']
        data_HI_fut_ssp245=HI_mod[modele][station]['ssp245']
        data_HI_fut_ssp585=HI_mod[modele][station]['ssp585']
    
        data_HI_hist_fut_ssp245=pd.concat([data_HI_hist,data_HI_fut_ssp245])
        data_HI_hist_fut_ssp585=pd.concat([data_HI_hist,data_HI_fut_ssp585])
        
        # annual maxima

        data_HI_ssp245=data_HI_hist_fut_ssp245.groupby(data_HI_hist_fut_ssp245.time.dt.year).max().HI.loc[period[0]:period[1]]
        data_HI_ssp585=data_HI_hist_fut_ssp585.groupby(data_HI_hist_fut_ssp585.time.dt.year).max().HI.loc[period[0]:period[1]]
    
        data_HI=np.concatenate((data_HI_ssp245.values.reshape(-1,1),data_HI_ssp585.values.reshape(-1,1)),axis=0)
            
        Y=data_HI
    
        #covariate
        
        data_covar_ssp245=df_covar[modele]['ssp245'].HI.loc[period[0]:period[1]]
        data_covar_ssp585=df_covar[modele]['ssp585'].HI.loc[period[0]:period[1]]
        data_covar=np.concatenate((data_covar_ssp245.values.reshape(-1,1),data_covar_ssp585.values.reshape(-1,1)),axis=0)
        
    
        X=data_covar
        
        # Averaged return levels as parameter
        
        RL =  np.array([RL_hist_ssp245[mod],RL_hist_ssp585[mod]])
        
        # GEV fit
        
        gevt=climextremes.fit_gev(Y,X, locationFun = 1, scaleFun = 1, returnValue = RL, 
                              optimArgs={"method" : "BFGS"},
                              initial={"location" : 40, "scale" : 1, "shape" : -0.01},
                              getParams=True, 
                              #bootSE=True
                             #getFit=True
                             )#likelihood maximisation : algo BFGS (Broyden-Fletcher-Goldfarb-Shanno)
        
        if gevt['info']['failure'][0]==0: # if convergence
            
            return_periods_ssp245=gevt['logReturnPeriod'][:(period[1]-period[0])+1,0]
            return_periods_ssp585=gevt['logReturnPeriod'][(period[1]-period[0])+1:,1]
            
            #storing RP for all models
                
            mat_ssp245[mod,:]=return_periods_ssp245
            mat_ssp585[mod,:]=return_periods_ssp585
                
        mod+=1
    
    # Saving results    
    
    gev_res['RP_mean']={}
    gev_res['RP_mean']['ssp245']=np.nanmean(np.exp(mat_ssp245), axis=0)
    gev_res['RP_mean']['ssp585']=np.nanmean(np.exp(mat_ssp585), axis=0)
    
    data_array_ssp245 = xr.Dataset(dict(
        
        RP =(["year"], gev_res['RP_mean']['ssp245']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RP_mean']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return periods of HI : year"),

        )
    
    data_array_ssp585 = xr.Dataset(dict(
        
        RP =(["year"], gev_res['RP_mean']['ssp585']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RP_mean']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return periods of HI : year"),

        )
    
    data_array_ssp245.to_netcdf(path=data_path+"GEV/{}y-RP_HI_1990-2100_ssp245_".format(RP)+station+'.nc')
    data_array_ssp585.to_netcdf(path=data_path+"GEV/{}y-RP_HI_1990-2100_ssp585_".format(RP)+station+'.nc')
    


# In[23]:


def calcul_conf_int_RL(station,period,RP):
    
    '''
    
    Compute multimodel confidence intervals of return levels
    The gev fit is done with both data from SSP2-4.5 and SSP5-8.5 scenarios,
    considering that gev parameters depend only on covariate values.
    
    In the paper: RP = 10 years, period=[1990,2100]
    
    This method returns:
    - Files storing bootstraped results in return levels, in both scenarios SSP2-4.5 and SSP5-8.5
    
    '''
        
    #multimodel results storage
    
    RL_results=[]
    ind_list=[]
    gev_res={}


    # GEV parameter storage
    
    param_gev=[]
    
    
    
    time_init=time.time()
    
    
    for modele_membre in Modeles :
        
        modele=modele_membre[:-9]
        
        #model results storage
        
        RL_results_mod=[]
        ind_list_mod=[]
        RL_results_hist_mod[modele]={}
        
        # modele data
        
        data_HI_hist=HI_mod[modele][station]['hist']
        data_HI_fut_ssp245=HI_mod[modele][station]['ssp245']
        data_HI_fut_ssp585=HI_mod[modele][station]['ssp585']
    
        data_HI_hist_fut_ssp245=pd.concat([data_HI_hist,data_HI_fut_ssp245])
        data_HI_hist_fut_ssp585=pd.concat([data_HI_hist,data_HI_fut_ssp585])
        
        #annual maxima
        
        data_HI_ssp245=data_HI_hist_fut_ssp245.groupby(data_HI_hist_fut_ssp245.time.dt.year).max().HI.loc[period[0]:period[1]]
        data_HI_ssp585=data_HI_hist_fut_ssp585.groupby(data_HI_hist_fut_ssp585.time.dt.year).max().HI.loc[period[0]:period[1]]
    
        data_HI=np.concatenate((data_HI_ssp245.values.reshape(-1,1),data_HI_ssp585.values.reshape(-1,1)),axis=0)
        
        Y=data_HI
    
        #covariate = HI
        data_covar_ssp245=df_covar[modele]['ssp245'].HI.loc[period[0]:period[1]]
        data_covar_ssp585=df_covar[modele]['ssp585'].HI.loc[period[0]:period[1]]
        data_covar=np.concatenate((data_covar_ssp245.values.reshape(-1,1),data_covar_ssp585.values.reshape(-1,1)),axis=0)
        
        #covariate = TX
        #data_covar=df_covar[modele][scenario].TX.loc[period[0]:period[1]]
    
        X=data_covar
        
        #bootstrap
        
        n=500 #number of iterations
    
        randind,randX,randY=bootstrap(X,Y,n)
        
        # GEV fits
        
        for i in range(n):
            
            gevt=climextremes.fit_gev(randY[i,:],randX[i,:], locationFun = 1, scaleFun = 1, returnPeriod = RP, 
                              optimArgs={"method" : "BFGS"},
                              initial={"location" : 40, "scale" : 1, "shape" : -0.01},
                              getParams=True, 
                              #bootSE=True
                             #getFit=True
                             )#likelihood maximisation : algo BFGS (Broyden-Fletcher-Goldfarb-Shanno)

    
            if gevt['info']['failure'][0]==0:# if convergence
            
                param_gev.append(gevt['mle'].tolist())
                
                #multimodel
                RL_results.append(gevt['returnValue'].tolist())
                ind_list.append(randind[i,:].tolist())
        
    
    #sort values in the good order
    
    RL_results=np.array(RL_results)
    ind_tab=np.array(ind_list)
    
    RL_res_sort=np.empty(RL_results.shape)
    RL_res_sort[:]=np.nan
    
    for i in range(np.shape(RL_res_sort)[0]):
        for j in range(np.shape(RL_res_sort)[1]):
            indij=ind_tab[i,j]
            RL_res_sort[i,indij]=RL_results[i,j]
            
    print('RL - number of series : ',RL_res_sort.shape[0])
    
    # saving results
    
    # samples to get uncertainties
    np.save(data_path+"RL_RP_HIX/{}y-RL_HI_1990-2100_ssp245_".format(RP)+station+'_all_samples.npy',RL_res_sort[:,:(period[1]-period[0])+1])
    np.save(data_path+"RL_RP_HIX/{}y-RL_HI_1990-2100_ssp585_".format(RP)+station+'_all_samples.npy',RL_res_sort[:,(period[1]-period[0])+1:])
    
    # quantile 0.05
    gev_res['RL_q5']={}
    gev_res['RL_q5']['ssp245']=np.nanquantile(RL_res_sort[:,:(period[1]-period[0])+1],0.05,axis=0)
    gev_res['RL_q5']['ssp585']=np.nanquantile(RL_res_sort[:,(period[1]-period[0])+1:],0.05,axis=0)
    
    #quantile 0.95
    gev_res['RL_q95']={}
    gev_res['RL_q95']['ssp245']=np.nanquantile(RL_res_sort[:,:(period[1]-period[0])+1],0.95,axis=0)
    gev_res['RL_q95']['ssp585']=np.nanquantile(RL_res_sort[:,(period[1]-period[0])+1:],0.95,axis=0)
    
    #GEV parameters
    np.save(data_path+"RL_RP_HIX/"+station+"_param_gev.npy",np.array(param_gev))
    
    data_array_ssp245_q5 = xr.Dataset(dict(
        
        RL =(["year"], gev_res['RL_q5']['ssp245']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RL_mean']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return levels of HI : °C"),

        )
    
    data_array_ssp585_q5 = xr.Dataset(dict(
        
        RL =(["year"], gev_res['RL_q5']['ssp585']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RL_q5']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return levels of HI : °C"),

        )
    
    
    data_array_ssp245_q95 = xr.Dataset(dict(
        
        RL =(["year"], gev_res['RL_q95']['ssp245']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RL_q5']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return levels of HI : °C"),

        )
    
    data_array_ssp585_q95 = xr.Dataset(dict(
        
        RL =(["year"], gev_res['RL_q95']['ssp585']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RL_q5']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return levels of HI : °C"),

        )
    
    data_array_ssp245_q5.to_netcdf(path=data_path+"GEV/{}y-RL_HI_1990-2100_ssp245_".format(RP)+station+'_q5.nc')
    data_array_ssp585_q5.to_netcdf(path=data_path+"GEV/{}y-RL_HI_1990-2100_ssp585_".format(RP)+station+'_q5.nc')
    
    
    data_array_ssp245_q95.to_netcdf(path=data_path+"GEV/{}y-RL_HI_1990-2100_ssp245_".format(RP)+station+'_q95.nc')
    data_array_ssp585_q95.to_netcdf(path=data_path+"GEV/{}y-RL_HI_1990-2100_ssp585_".format(RP)+station+'_q95.nc')
    
    
    time_end=time.time()
    
    print('RL - Time for {} iterations : {:.3f} s'.format(n,time_end-time_init))


# In[24]:


def calcul_conf_int_RP(station,period,RL_hist_ssp245,RL_hist_ssp585,RP):
    
    '''
    
    Compute multimodel confidence intervals of return levels
    The gev fit is done with both data from SSP2-4.5 and SSP5-8.5 scenarios,
    considering that gev parameters depend only on covariate values.
    
    In the paper: RP = 10 years, period=[1990,2100], RL_hist_ssp245
    and RL_hist_ssp585 are averaged historical return levels
    
    This method returns:
    - Files storing bootstraped results in return periods, in both scenarios SSP2-4.5 and SSP5-8.5
    
    '''
    
    
    #multimodel results storage
    
    RP_results=[]
    ind_list=[]
    
    time_init=time.time()
    
    mod=0
    
    for modele_membre in Modeles :
                
        modele=modele_membre[:-9]
        
        #model data

        data_HI_hist=HI_mod[modele][station]['hist']
        data_HI_fut_ssp245=HI_mod[modele][station]['ssp245']
        data_HI_fut_ssp585=HI_mod[modele][station]['ssp585']
    
        data_HI_hist_fut_ssp245=pd.concat([data_HI_hist,data_HI_fut_ssp245])
        data_HI_hist_fut_ssp585=pd.concat([data_HI_hist,data_HI_fut_ssp585])
        
        #annual maxima
        
        data_HI_ssp245=data_HI_hist_fut_ssp245.groupby(data_HI_hist_fut_ssp245.time.dt.year).max().HI.loc[period[0]:period[1]]
        data_HI_ssp585=data_HI_hist_fut_ssp585.groupby(data_HI_hist_fut_ssp585.time.dt.year).max().HI.loc[period[0]:period[1]]
    
        data_HI=np.concatenate((data_HI_ssp245.values.reshape(-1,1),data_HI_ssp585.values.reshape(-1,1)),axis=0)
        
        Y=data_HI
    
        #covariate
        
        data_covar_ssp245=df_covar[modele]['ssp245'].HI.loc[period[0]:period[1]]
        data_covar_ssp585=df_covar[modele]['ssp585'].HI.loc[period[0]:period[1]]
        data_covar=np.concatenate((data_covar_ssp245.values.reshape(-1,1),data_covar_ssp585.values.reshape(-1,1)),axis=0)
        
    
        X=data_covar
        
        #bootstrap
        
        n=500 # number of iterations
    
        randind,randX,randY=bootstrap(X,Y,n)
        
        RL =  np.array([RL_hist_ssp245[mod],RL_hist_ssp585[mod]])
        
        # GEV fits
        
        for i in range(n):
            
            gevt=climextremes.fit_gev(randY[i,:],randX[i,:], locationFun = 1, scaleFun = 1, 
                                returnValue = RL, 
                              optimArgs={"method" : "BFGS"},
                              initial={"location" : 40, "scale" : 1, "shape" : -0.01},
                              getParams=True, 
                              #bootSE=True
                             #getFit=True
                             )#likelihood maximisation : algo BFGS (Broyden-Fletcher-Goldfarb-Shanno)

    
            if gevt['info']['failure'][0]==0:# if convergence
                
                
                RP_results.append(np.exp(gevt['logReturnPeriod']).tolist())
                ind_list.append(randind[i,:].tolist())
                
        mod+=1
        
    #sort values in the good order
    
    RP_results=np.array(RP_results)
    ind_tab=np.array(ind_list)    
            
    RP_res_sort=np.empty(RP_results.shape)
    RP_res_sort[:]=np.nan
    
    for i in range(np.shape(RP_res_sort)[0]):
        for j in range(np.shape(RP_res_sort)[1]):
            indij=ind_tab[i,j]
            RP_res_sort[i,indij]=RP_results[i,j]
    
    print('RP - number of series : ',RP_res_sort.shape[0])
    
    # saving results
    
    gev_res={}
        
    # samples, quantile 0.05 and quantile 0.95
    
    np.save(data_path+"RL_RP_HIX/{}y-RP_HI_1990-2100_ssp245_".format(RP)+station+'_all_samples.npy',RP_res_sort[:,:(period[1]-period[0])+1,0])
    np.save(data_path+"RL_RP_HIX/{}y-RP_HI_1990-2100_ssp585_".format(RP)+station+'_all_samples.npy',RP_res_sort[:,(period[1]-period[0])+1:,1])
    
    gev_res['RP_q5']={}
    gev_res['RP_q5']['ssp245']=np.nanquantile(RP_res_sort[:,:(period[1]-period[0])+1,0],0.05,axis=0)
    gev_res['RP_q5']['ssp585']=np.nanquantile(RP_res_sort[:,(period[1]-period[0])+1:,1],0.05,axis=0)
    
    gev_res['RP_q95']={}
    gev_res['RP_q95']['ssp245']=np.nanquantile(RP_res_sort[:,:(period[1]-period[0])+1,0],0.95,axis=0)
    gev_res['RP_q95']['ssp585']=np.nanquantile(RP_res_sort[:,(period[1]-period[0])+1:,1],0.95,axis=0)
    
    data_array_ssp245_q5 = xr.Dataset(dict(
        
        RP =(["year"], gev_res['RP_q5']['ssp245']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RP_q5']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return periods of HI : year"),

        )
    
    data_array_ssp585_q5 = xr.Dataset(dict(
        
        RP =(["year"], gev_res['RP_q5']['ssp585']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RP_q5']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return periods of HI : year"),

        )
    
    
    data_array_ssp245_q95 = xr.Dataset(dict(
        
        RP =(["year"], gev_res['RP_q95']['ssp245']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RP_q5']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return periods of HI : year"),

        )
    
    data_array_ssp585_q95 = xr.Dataset(dict(
        
        RP =(["year"], gev_res['RP_q95']['ssp585']),

        ),
                               
        coords=dict(

        year=("year", data_covar_ssp245.index.values[:np.shape(gev_res['RP_q5']['ssp245'])[0]]),
        
        ),
        
        attrs=dict(units="Return periods of HI : year"),

        )
    
    data_array_ssp245_q5.to_netcdf(path=data_path+"GEV/{}y-RP_HI_1990-2100_ssp245_".format(RP)+station+'_q5.nc')
    data_array_ssp585_q5.to_netcdf(path=data_path+"GEV/{}y-RP_HI_1990-2100_ssp585_".format(RP)+station+'_q5.nc')

    data_array_ssp245_q95.to_netcdf(path=data_path+"GEV/{}y-RP_HI_1990-2100_ssp245_".format(RP)+station+'_q95.nc')
    data_array_ssp585_q95.to_netcdf(path=data_path+"GEV/{}y-RP_HI_1990-2100_ssp585_".format(RP)+station+'_q95.nc')
    
    
    time_end=time.time()
    
    print('RP - Time for {} iterations : {:.3f} s'.format(n,time_end-time_init))
    


# ## RUN SCRIPTS

# In[25]:


for station in Stations:
    for RP in [10]:
        calcul_extreme_stats_RL(station,[1990,2100],RP)
        calcul_conf_int_RL(station,[1990,2100],RP)


# In[26]:


for station in Stations:
    for RP in [10]:
        RL_hist_ssp245=xr.open_dataset(data_path+'GEV/'+str(RP)+'y-RL_HI_hist_ssp245_'+station+'.nc')
        RL_hist_ssp585=xr.open_dataset(data_path+'GEV/'+str(RP)+'y-RL_HI_hist_ssp585_'+station+'.nc')
    
        calcul_extreme_stats_RP(station,[1990,2100],RL_hist_ssp245.RL.values,RL_hist_ssp585.RL.values,RP)
        calcul_conf_int_RP(station,[1990,2100],RL_hist_ssp245.RL.values,RL_hist_ssp585.RL.values,RP)


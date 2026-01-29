#! /usr/bin/python
# =============================================================================

# Content : CDF-t correction for minimal daily temperature in future period (2015-2100)
# Modeles : CMIP6

# =============================================================================
# -- Python packages
# -----------------------------------------------------------------------------
from   cdo import *
import datetime
import numpy as np
import xarray as xr
import numpy.ma as ma
import time as tm 
import os, stat, glob
from   netCDF4 import Dataset
from   statsmodels.distributions.empirical_distribution import ECDF
from joblib import Parallel, delayed
cdo = Cdo()

# =============================================================================
#  -- Main variables
# -----------------------------------------------------------------------------
list_scenarios= ["ssp245","ssp585"]

list_models   = ["ACCESS-CM2","CanESM5","CNRM-ESM2-1","CNRM-CM6-1","CNRM-CM6-1-HR","EC-Earth3","EC-Earth3-CC","EC-Earth3-Veg-LR",
                 "IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR","NorESM2-MM"]

list_vars    = ["tasmin"]

list_stations=["Gillot","Pamandzi","Tromelin"]

unit         = 'dgC'
BCM          = 'BC_3M'# 3 sliding months
cdf_thres_lower = .001 # to avoid outliers when correcting
cdf_thres_upper = .999 # to avoid outliers when correcting
HIST_period  =  "1985-2014" # Historical period
OBS_period = "1985-2014" # Reference period
HST_ref      = "1985/2014" # common years between Historical and reference
MOD_vartime  = "time"

###########################################################################

def corr_over_models(modele):
    
            variable = "_".join((MOD_var,"day"))
            modele_alias = modele
        
            
            FUTUR_period = "2015-2100"
            FUT_ref      = '2015/2100' 
            FUT_ref_stc  = '2071/2100' # Static future period
            HST_FUT_period = '2000-2100' # when combining historical and future
        
        

            if (modele=="MIROC-ES2L"):
              version      = "r1i1p1f2_gn"
            elif (modele=="ACCESS-CM2"):
              version      = "r2i1p1f1_gn"
            elif (modele=="CNRM-CM6-1" or modele=="CNRM-CM6-1-HR" or modele=="CNRM-ESM2-1"):
              version      = "r1i1p1f2_gr"
            elif (modele=="IPSL-CM6A-LR" or modele=="EC-Earth3" or modele=="EC-Earth3-Veg-LR" or modele=="EC-Earth3-CC" or modele=="KACE-1-0-G"):
                version      = "r1i1p1f1_gr"
            else:
                version    = "r1i1p1f1_gn"
              
            if (modele=="IPSL-CM6A-LR" or modele=="KACE-1-0-G"):
              date_ref = 'days since 2015-01-01 00:00:00'
            elif (modele=="NorESM2-MM"):
              date_ref = 'days since 1-01-01 00:00:00'
            else:
              date_ref = 'days since 1850-01-01 00:00:00'
              
            if (modele=="CanESM5" or modele=="NorESM2-MM"):
                 calendar_att='365-day'
            elif (modele=="KACE-1-0-G"):
                calendar_att='360-day'
            else :
                calendar_att='standard'
        
            # =============================================================================
            #for ind_scen in range(0,len(list_scenarios)):
            scenario = list_scenarios[ind_scen]
          
            HIST         = "_".join((modele,"historical",version[:-3])) 
            FUTUR        = "_".join((modele,scenario,version[:-3]))  
            FUTUR_BC     = "_".join((FUTUR,"BC_3M"))
            
            # =============================================================================
            #  -- Folders
            # -----------------------------------------------------------------------------
            
            dirin        = "./input_cdft/"+MOD_var+"/" # folder for raw data
            dirout       = "./output_cdft/"+MOD_var+"/" # folder for corrected data          
            tempo_dir    = './output_cdft/tmp/'  # temporary folder
            OBS_dirin    = "./Obs/" 
            HIST_dirin   = dirin
            FUTUR_dirin  = dirin
            FUTUR_dirout = dirout
            cdo = Cdo(tempdir = tempo_dir)
            
            # =============================================================================
            # -- Files
            # ----------------------------------------------------------------------------- 
            prefin_OBS   = "_".join((station,'obs', OBS_period, MOD_var))
            prefin_HIST   = "_".join((station,'hist',MOD_var,modele))
            prefin_FUTUR  = "_".join((station,'futur',scenario,MOD_var,modele))
            
            prefout       = "_".join((station,FUTUR_BC,variable))
            prefout_FUTUR = "_".join((station,modele_alias,BCM,variable,scenario,FUTUR_period))  
            ficoutf_FUTUR = os.path.join(FUTUR_dirout,prefout_FUTUR) + '.nc'
            print('prefout_FUTUR : ' + prefout_FUTUR)
            
            # =============================================================================
            # -- Read netCDF HIST and merge daily data input file
            # -----------------------------------------------------------------------------
            start0 = tm.time()
            print ('\n*******************************************************************')
            print ('Start: ', datetime.datetime.fromtimestamp(tm.time()).strftime('%c'))
            print ('Inputs files reading ......')   
            obs_day = os.path.join(OBS_dirin,prefin_OBS) + '.nc' # reference
            hst_day = os.path.join(HIST_dirin,prefin_HIST) + '.nc' # modele historical
            fut_day = os.path.join(FUTUR_dirin,prefin_FUTUR) + '.nc' # modele future
            print('Merged  input file is reading from ' + obs_day)  
            print('History input file is reading from ' + hst_day)  
            print('Future  input file is reading from ' + fut_day)

            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            
            # =============================================================================
            # Merge history (1985-2014) and scenario (2015-2100) data
            # -----------------------------------------------------------------------------
            new_fut_day = FUTUR_dirin + '/' + MOD_var + '_' +station+'_'+scenario + '_'+modele+'_'+ HST_FUT_period + '.nc'
            print (MOD_var, scenario, 'History_', HST_ref, ' and Scenario_', FUTUR_period,' and Modele_',modele, ' merging ......')
            start1 = tm.time()
            cdo.mergetime(input = cdo.selyear(HST_ref, input = hst_day) + ' ' + 
                                  cdo.selyear(FUT_ref, input = fut_day), 
                                  output = new_fut_day)
            
            end1 = tm.time()
            print (MOD_var, scenario, 'History_', HST_ref, ' and Scenario_', FUTUR_period, ' merging is finished,'
                  '\nDuration:', (end1 - start1)/60, 'minutes', sep=' ') 
            print ('-------------------------------------------------------------------\n')
            
            # -----------------------------------------------------------------------------
            # Moving year range
           
            window_inf  = list(np.arange(2000,2086,1)) # 15 years before adjusted year list
            window_cnt  = list(np.arange(2015,2101,1)) # adjusted year list
            window_sup  = list(np.arange(2029,2115,1)) # 15 years after adjusted year list
            
            
            # =============================================================================
            # -- Initialisation
            # -----------------------------------------------------------------------------
            for i in range(len(window_cnt)):
               start2   = tm.time()
               start_yr = str(window_inf[i])
               corr_yr  = str(window_cnt[i]) # adjusted year
               end_yr   = str(window_sup[i])
            
               moi   = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
               N1MON = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 
            	    'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
            
               FUT_ref_mvg     = '/'.join((start_yr,end_yr))
               prefout_corr_yr = "_".join((prefout,corr_yr))
               
               print (scenario, corr_yr, BCM, MOD_var,'Bias Correction ...')
               print ('------------------------------------------\n')
            
               for counter in range(len(moi)):
                  start3 = tm.time()
                  print ('FUT ', scenario, corr_yr, N1MON[counter%12], 
                                     'and her window data selection ...')
                  
                  # Select future monthly data for adjusted year 'corr_yr'
                  fut_mon = cdo.selmon(moi[counter%12], 
                                       input   = '-selyear,' + corr_yr + ' %s' % (fut_day), 
                                       options = '-f nc', returnCdf = True)

                  
                  # Select future window for adjusted year 'corr_yr'
                  if (corr_yr <= '2085'):
                     fut_window = cdo.selmon(moi[(counter-1)%12], moi[counter%12], 
                                             moi[(counter+1)%12], 
                                             input   = '-selyear,' + FUT_ref_mvg + ' %s' 
                                                         % (new_fut_day), 
                                             options = '-f nc', returnCdf = True) # 3 months over 30 years
                  else: #after 2085 we need the static period 2071-2100
                     fut_window = cdo.selmon(moi[(counter-1)%12], moi[counter%12], 
                                             moi[(counter+1)%12], 
                                             input   = '-selyear,' + FUT_ref_stc + ' %s' 
                                                         % (fut_day), 
                                             options = '-f nc', returnCdf = True)
                     
                  # Select history window for reference period (1985-2014)
                  hst_window = cdo.selmon(moi[(counter-1)%12], moi[counter%12], 
                                          moi[(counter+1)%12], 
                                          input = '-selyear,' + HST_ref + ' %s' % (hst_day), 
                                          options = '-f nc', returnCdf = True)
                  
                  # Select obs window for reference period : 1985/2014
                  mrg_window = cdo.selmon(moi[(counter-1)%12], moi[counter%12], 
                                          moi[(counter+1)%12], 
                                          input = '-selyear,' + HST_ref + ' %s' % (obs_day), 
                                          options = '-f nc', returnCdf = True)
            
                  # =======================================================================
                  # -- Output can be read directory as a netcdf file
                  # -----------------------------------------------------------------------
                  times = fut_mon.variables[MOD_vartime][:] 
               
                  tas_fut_mon    = fut_mon.variables[MOD_var][:]
                  tas_fut_window = fut_window.variables[MOD_var][:]
                  tas_hst_window = hst_window.variables[MOD_var][:]
                  tas_obs_window = mrg_window.variables[MOD_var][:]
               
                  # =======================================================================
                  # -- Create new Array to initialize BC_pr_fut_mon 
                  # -----------------------------------------------------------------------
                  BC_tas_fut_mon = np.float32(-99+np.zeros((len(times)))) 
                  
                  # =======================================================================
                  #                       MAIN PROGRAM
                  # =======================================================================
                  # -- Loops for lats, lons and times variations
                  # -----------------------------------------------------------------------
                  print ('  ', scenario, corr_yr, N1MON[counter%12], MOD_var, 
            					'Bias Correction ...')
               
                  ecdf = ECDF(tas_fut_window[:,0,0])
                  fut_mon_ecdf = ecdf(tas_fut_mon[:,0,0])
                  fut_mon_ecdf[fut_mon_ecdf > cdf_thres_upper] = cdf_thres_upper
                  fut_mon_ecdf[fut_mon_ecdf < cdf_thres_lower] = cdf_thres_lower
                  
                  # Inverse CDF estimated from observed and historical values
                  obs_mon_quantil = np.nanquantile(tas_obs_window[:], fut_mon_ecdf)
                  hst_mon_quantil = np.nanquantile(tas_hst_window[:,0,0], fut_mon_ecdf)
                  
                  # Correction coefficient
                  corr_tas_quant = obs_mon_quantil - hst_mon_quantil
                  
                  # Correction
                  BC_tas_fut_mon[:] = tas_fut_mon[:,0,0] + corr_tas_quant
                  
                  end3 = tm.time()
                  print('  ', scenario, corr_yr, N1MON[counter%12], MOD_var,
            		'Bias Correction is finished,' 
                        'duration:', (end3 - start3)/60, 'minutes', sep=' ') 
               
                  # =======================================================================
                  # -- Save Dataset to NetCDF
                  # -----------------------------------------------------------------------
                  prefout_FUT   = "_".join((prefout_corr_yr, N1MON[counter%12]))
                  ficout_FUTUR  = os.path.join(tempo_dir,prefout_FUT) + '_g.nc'
                  ficoutg_FUTUR = os.path.join(tempo_dir,prefout_FUT) + '.nc'
                  print ('   Now creating ' + ficout_FUTUR)
            
                  # Create writable netcdf file
                  dataset = Dataset(ficout_FUTUR, 'w', format='NETCDF4_CLASSIC') 
               
                  #Define a set of dimensions used for our variables
                  time = dataset.createDimension('time', len(times))
               
                  # Create coordinate variables for 3-dimensions
                  #Dataset.createVariable(<var_id>, <type>, <dimensions>)
                  time_counter = dataset.createVariable('time', np.float64, ('time',))
                  tas_corr      = dataset.createVariable(MOD_var + '_corr', np.float32, 
                                                        ('time',))
               
                  # Global Attributes
                  dataset.description = 'BIAS CORRECTION for city' + station
                  dataset.history     = 'Created ' + tm.ctime(tm.time())
               
                  # Variable Attributes
                  tas_corr.units         = unit
                  tas_corr.long_name     = long_name
                  tas_corr.fill_value    = -99
                  tas_corr.missing_value = -99
                  time_counter.units    = date_ref 
                  time_counter.calendar = calendar_att
                  time_counter.axis     = 'T'
               
                  # Store data into the netcdf file
                  time_counter[:] = times
                  tas_corr[:]  = xr.DataArray(BC_tas_fut_mon, coords=[times], 
                                                 dims=['time'])  
                  dataset.close()   # and the file is now written!
                  
                  print('   --------------------------------------------------------\n')
                              
               # ==========================================================================
               # Merging all monthly files
               # --------------------------------------------------------------------------
               ficout_corr_yr = os.path.join(tempo_dir,prefout_corr_yr)
               cdo.mergetime(input = tempo_dir + '/' + prefout_corr_yr + '_*.nc', 
                             output = ficout_corr_yr + '.nc')
               
               end2 = tm.time()
               print (scenario, corr_yr, BCM, MOD_var, 'Bias correction is finished,', 
            		'duration:', (end2 - start2)/3600, 'hours') 
               print ('Now aggregating the 12 files into ' + prefout_corr_yr + '.nc')
               print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            
               # Remove all files in tempo_dir directory
               test = tempo_dir + '/' + prefout_corr_yr + '_*.nc'
               r = glob.glob(test)
               for i in r:
                   os.remove(i)
            
            end0 = tm.time()
            print (scenario, BCM, MOD_var, 'Bias Correction is finished successfuly\n'
                   'Duration:' , (end0 - start0)/3600, 'hours', sep=' ')
            print ('End: ', datetime.datetime.fromtimestamp(tm.time()).strftime('%c'))
        
            # ==========================================================================
            # Merging all monthly files
            # --------------------------------------------------------------------------
            print ('Now aggregating the yearly files into ' + ficoutf_FUTUR)
            print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            cdo.mergetime(input = tempo_dir+ '/' + prefout + '*.nc', 
                             output = ficoutf_FUTUR)
            # Remove all files in tempo_dir directory
            test = tempo_dir + '/' +  prefout + '*.nc'
            r = glob.glob(test)
            for i in r:
                os.remove(i)
            print ('*******************************************************************\n') 
    

###########################################################################

# Job's parallelisation over models

MOD_var = list_vars[0]
long_name = '_'.join(('Daily', MOD_var, 'after Bias Correction'))

for ind_scen in range(len(list_scenarios)):
    
    for station in list_stations :
    
        Parallel(n_jobs=len(list_models))(delayed(corr_over_models)(modele) for modele in list_models)       

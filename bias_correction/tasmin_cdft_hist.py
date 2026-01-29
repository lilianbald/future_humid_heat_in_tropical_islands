#! /usr/bin/python
# =============================================================================

# Content : CDF-t correction for minimal daily temperature in Historical period (1985-2014)
# Modeles : CMIP6

# =============================================================================
# -- Python packages
# -----------------------------------------------------------------------------
from   cdo import *
import datetime
import numpy as np
import xarray as xr
import time as tm
import os, glob
from   netCDF4 import Dataset
from   statsmodels.distributions.empirical_distribution import ECDF
from joblib import Parallel, delayed
cdo = Cdo()

# =============================================================================
#  -- Main variables
# -----------------------------------------------------------------------------
list_models   = ["ACCESS-CM2","CanESM5","CNRM-ESM2-1","CNRM-CM6-1","CNRM-CM6-1-HR","EC-Earth3","EC-Earth3-CC","EC-Earth3-Veg-LR",
                 "IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR","NorESM2-MM"]

list_vars   = ["tasmin"]

list_stations=["Gillot","Pamandzi","Tromelin"]

unit         = 'dgC'
BCM          = 'BC_3M'# 3 sliding months
cdf_thres_lower = .001 # to avoid outliers when correcting
cdf_thres_upper = .999 # to avoid outliers when correcting
HIST_period  =  "1985-2014" # Historical period
OBS_period = "1985-2014" # Reference period
FUTUR_period = "2015-2100"
HST_ref      = "1985/2014" # Common years between Historical and reference
HST_ref_stc  = "1985/2014" # Static period
MOD_vartime  = "time" 

##########################################################################
def corr_over_models(modele):
        modele_alias = modele
        # =============================================================================
        # -- Make our variables change here
        # ----------------------------------------------------------------------------
        scenario  = 'ssp245'
        
        # =============================================================================
        
        variable = "_".join((MOD_var,"day"))
        

        if (modele=="MIROC-ES2L"):
          version      = "r1i1p1f2_gn"
       
        elif (modele=="ACCESS-CM2):
          version      = "r2i1p1f1_gn"
        elif (modele=="CNRM-CM6-1" or modele=="CNRM-CM6-1-HR" or modele=="CNRM-ESM2-1"): 
        elif (modele=="IPSL-CM6A-LR" or modele=="EC-Earth3" or modele=="EC-Earth3-Veg-LR" or modele=="EC-Earth3-CC" or modele=="KACE-1-0-G"):
            version      = "r1i1p1f1_gr"
        else:
            version    = "r1i1p1f1_gn"
            
        if ( modele=="NorESM2-MM"):
          date_ref = 'days since 1-01-01 00:00:00'
        else:  
          date_ref = 'days since 1850-01-01 00:00:00'
          
        if (modele=="CanESM5" or modele=="NorESM2-MM"):
            calendar_att='365-day'
        elif (modele=="KACE-1-0-G"):
            calendar_att='360-day'
        else :
            calendar_att='standard'
    
        HIST         = "_".join((modele,"historical",version[:-3])) 
        HIST_BC      = "_".join((HIST,"BC_3M"))
        FUTUR        = "_".join((modele,scenario,version[:-3])) 
      
        # =============================================================================
        #  -- Folders
        # -----------------------------------------------------------------------------
        dirin        = "./input_cdft/"+MOD_var+"/" #folder for raw data
        dirout       = "./output_cdft/"+MOD_var+"/"#folder for corrected data
        tempo_dir    = './output_cdft/tmp/'  # temporary folder
        OBS_dirin    = "./Obs/" 
        HIST_dirin   = dirin
        HIST_dirout  = dirout
        FUTUR_dirin  = dirin # where we can find also data over 2015-2100
        cdo = Cdo(tempdir = tempo_dir)
        
        # =============================================================================
        # -- Files
        # ----------------------------------------------------------------------------- 
        prefin_OBS  = "_".join((station,'obs', OBS_period, MOD_var))
        prefin_HIST   = "_".join((station,'hist',MOD_var,modele))
        prefin_FUTUR  = "_".join((station,'futur',scenario,MOD_var,modele))
        
        prefout      = "_".join((station,HIST_BC,variable))
        prefout_HIST = "_".join((station,modele_alias,BCM,variable,HIST_period))  
        ficoutf_HIST = os.path.join(HIST_dirout,prefout_HIST) + '.nc'
        print('prefout_HIST : ' + prefout_HIST)
      
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
        print('Obs input file is reading from ' + obs_day)  
        print('History input file is reading from ' + hst_day)  
        print('Future  input file is reading from ' + fut_day)
        
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        
        # =============================================================================
        # Merge history (2000-2014) and scenario (2015-2100) data
        # -----------------------------------------------------------------------------
        new_hst_day = FUTUR_dirin + '/' + MOD_var +'_'+station+'_'+modele+ '_HIST_1985-2030.nc'
        print (MOD_var, scenario, 'History_', HIST_period,modele, 'and Scenario_2015-2030 merging ......')
        start1 = tm.time()
        cdo.mergetime(input = cdo.selyear(HST_ref, input = hst_day) + ' ' + 
                              cdo.selyear('2015/2030', input = fut_day), 
                              output = new_hst_day)
        
        end1 = tm.time()
        print (MOD_var,'History_' , HIST_period, 'and Scenario_2015-2030 merging is finished,'
              '\nDuration:', (end1 - start1)/60, 'minutes', sep=' ') 
        print ('-------------------------------------------------------------------\n')
        
        # -----------------------------------------------------------------------------
        # Moving year range
        window_inf  = list(np.arange(1970,2000,1))# 15 years before adjusted year list
        window_cnt  = list(np.arange(1985,2015,1))# adjusted year list
        window_sup  = list(np.arange(1999,2029,1))# 15 years after adjusted year list
        
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
        
           HST_ref_mvg     = '/'.join((start_yr,end_yr))
           prefout_corr_yr = "_".join((prefout,corr_yr)) 
           
           print ('HIST ', corr_yr, BCM, MOD_var,'Bias Correction ...')
           print ('------------------------------------------\n')
        
           for counter in range(len(moi)):
              start3 = tm.time()
              print ('HIST ', corr_yr, N1MON[counter%12], 
                                 'and her window data selection ...')
              
              # Select future monthly data for adjusted year 'corr_yr'
              hst_mon = cdo.selmon(moi[counter%12], 
                                   input   = '-selyear,' + corr_yr + ' %s' % (hst_day), 
                                   options = '-f nc', returnCdf = True)# mois d'une année
              
              # Select future window for adjusted year 'corr_yr'
              if (corr_yr >= '2001'):# we need to include future period (>2015) in the sliding period
                 hst_window = cdo.selmon(moi[(counter-1)%12], moi[counter%12], 
                                         moi[(counter+1)%12], 
                                         input   = '-selyear,' + HST_ref_mvg + ' %s' 
                                                     % (new_hst_day), 
                                         options = '-f nc', returnCdf = True)# 3 months over 30 years
              else: #before 2001 we don't need future period, just use 1985-2014
                 hst_window = cdo.selmon(moi[(counter-1)%12], moi[counter%12], 
                                         moi[(counter+1)%12], 
                                         input   = '-selyear,' + HST_ref_stc + ' %s' 
                                                     % (new_hst_day), 
                                         options = '-f nc', returnCdf = True)
              
              # Select history window for reference period : 1985/2014
              hst_window_ref = cdo.selmon(moi[(counter-1)%12], moi[counter%12], 
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
              times = hst_mon.variables[MOD_vartime][:] 
           
              tas_hst_mon    = hst_mon.variables[MOD_var][:]
              tas_hst_window = hst_window.variables[MOD_var][:] #3 months
              tas_hst_window_ref = hst_window_ref.variables[MOD_var][:]#1985-2014
              tas_obs_window = mrg_window.variables[MOD_var][:]
           
              # =======================================================================
              # -- Create new Array to initialize BC_pr_fut_mon 
              # -----------------------------------------------------------------------
              BC_tas_hst_mon = np.float32(-99+np.zeros((len(times)))) 
              
              # =======================================================================
              #                       MAIN PROGRAM
              # =======================================================================
              # -- Loops for lats, lons and times variations
              # -----------------------------------------------------------------------
              print (' HIST ', corr_yr, N1MON[counter%12], MOD_var, 
        					'Bias Correction ...')
              # Nonexceedance probability associated with the value at time t
              ecdf = ECDF(tas_hst_window[:,0,0])
              hst_mon_ecdf = ecdf(tas_hst_mon[:,0,0])
              hst_mon_ecdf[hst_mon_ecdf > cdf_thres_upper] = cdf_thres_upper # we don't exceed cdf_thres_upper
              hst_mon_ecdf[hst_mon_ecdf < cdf_thres_lower] = cdf_thres_lower # we don't take values below cdf_thres_lower
                       
              # Inverse CDF estimated from observed and historical values
              obs_mon_quantil = np.nanquantile(tas_obs_window[:], hst_mon_ecdf)
              hst_mon_quantil = np.nanquantile(tas_hst_window_ref[:,0,0], hst_mon_ecdf)
                       
              # Correction coefficient
              corr_tas_quant = obs_mon_quantil - hst_mon_quantil
                       
              # Correction
              BC_tas_hst_mon[:] = tas_hst_mon[:,0,0] + corr_tas_quant
              
              end3 = tm.time()
              print(' HIST ', corr_yr, N1MON[counter%12], MOD_var,
        		'Bias Correction is finished,' 
                    'duration:', (end3 - start3)/60, 'minutes', sep=' ') 
           
              # =======================================================================
              # -- Save Dataset to NetCDF
              # -----------------------------------------------------------------------
              prefout_HST   = "_".join((prefout_corr_yr, N1MON[counter%12]))
              ficout_HIST  = os.path.join(tempo_dir,prefout_HST) + '_g.nc'
              ficoutg_HIST = os.path.join(tempo_dir,prefout_HST) + '.nc'
              print ('   Now creating ' + ficout_HIST)
        
              # Create writable netcdf file
              dataset = Dataset(ficout_HIST, 'w', format='NETCDF4_CLASSIC') 
           
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
              tas_corr[:]  = xr.DataArray(BC_tas_hst_mon, coords=[times], 
                                             dims=['time'])  
              dataset.close()   # and the file is now written!
              
              print('   --------------------------------------------------------\n')
                      
           # ==========================================================================
           # Merging all monthly files
           # --------------------------------------------------------------------------
           ficout_corr_yr = os.path.join(tempo_dir,prefout_corr_yr)
           cdo.mergetime(input = tempo_dir + prefout_corr_yr + '_*.nc', 
                         output = ficout_corr_yr + '.nc')
           
           end2 = tm.time()
           print ('HIST ', corr_yr, BCM, MOD_var, 'Bias correction is finished,', 
        		'duration:', (end2 - start2)/3600, 'hours') 
           print ('Now aggregating the 12 files into ' + prefout_corr_yr + '.nc')
           print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        
           # Remove all files in tempo_dir directory
           test = tempo_dir + '/' + prefout_corr_yr + '_*.nc'
           r = glob.glob(test)
           for i in r:
               os.remove(i)
        
        end0 = tm.time()
        print ('HIST ', BCM, MOD_var, 'Bias Correction is finished successfuly\n'
               'Duration:' , (end0 - start0)/3600, 'hours', sep=' ')
        print ('End: ', datetime.datetime.fromtimestamp(tm.time()).strftime('%c'))
    
        # ==========================================================================
        # Merging all monthly files
        # --------------------------------------------------------------------------
        print ('Now aggregating the yearly files into ' + ficoutf_HIST)
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        cdo.mergetime(input = tempo_dir+ '/' + prefout + '*.nc', 
                         output = ficoutf_HIST)
        # Remove all files in tempo_dir directory
        test = tempo_dir + '/' +  prefout + '*.nc'
        r = glob.glob(test)
        for i in r:
            os.remove(i)
        print ('*******************************************************************\n') 

# ##########################################################################

# Job's parallelisation over models

for ind_var in range(0,len(list_vars)):
  MOD_var      = list_vars[ind_var]
  long_name    = '_'.join(('Daily', MOD_var, 'after Bias Correction'))
  
  for station in list_stations:
      
      Parallel(n_jobs=len(list_models))(delayed(corr_over_models)(modele) for modele in list_models)

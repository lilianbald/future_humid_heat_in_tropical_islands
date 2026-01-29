# future_humid_heat_in_tropical_islands
This repository contains the data and codes whose associated results are summarised in the article "Small tropical islands also exposed to extreme humid heat by the end of the century", Bald et al.

## Title
Small tropical islands also exposed to extreme humid heat by the end of the century

## Authors
Lilian Bald*, Ali Belmadani, Marie-Dominique Leroux, Olivier Pannekoucke, Agathe Gentric, Saïd Qasmi. 

*: lilian.bald@meteo.fr

## Abstract
Extreme humid heat is projected to increase markedly by the end of the century, particularly in the tropics. Yet, the case of small islands is ambiguous because the available studies are based on coarse global climate models (GCMs) that ignore them. Here, statistically downscaled, bias-corrected GCMs at weather station sites are used to assess present and future extreme humid heat conditions on various tropical islands. The latter are estimated with the Heat Index (HI) combining temperature and humidity. Similarly to continents, extreme humid heat intensity is projected to reach particularly dangerous levels by the end of the century, with higher HI values for islands located closer to the equator. Longer, more common humid heatwaves are also notably projected to increase in frequency, with more pronounced increases equatorward, because of the lower seasonal variability of HI. Such severe conditions threaten the lives of millions of island inhabitants across the tropics.

## Data

Data is available on Zenodo: (DOI).

Three folders:

- ALADIN: contains bias-corrected data from RCM ALADIN in station Gillot.
  - ALADIN/GEV: contains the 10-year HIX return levels for Gillot, deduced from ALADIN.
- CMIP6: contains bias-corrected data from CMIP6 models.
  - CMIP6/covar: contains the covariate values for each CMIP6 model, each scenario and for each domain (North Atlantic -- "_Guyane.csv", Southwest Indian -- "_Reunion.csv", South Pacific -- "_NC.csv").
  - CMIP6/GEV: 10-year return levels for each model, each station, each scenario. "all_samples" are needed to compute confidence intervals.
  - CMIP6/heatwaves: characteristics (intensity, duration, frequency) of HHW1 and HHW2 for each station and scenario, in periods 1990-2000 and 2090-2100.
  - CMIP6/HIM: daily mean HI CMIP6 data for each station, each model, each scenario, over periods 1990-2014 and 2015-2100.
  - CMIP6/HIN: daily minimum HI CMIP6 data for each station, each model, each scenario, over periods 1990-2014 and 2015-2100.
  - CMIP6/HIX: daily maximum HI CMIP6 data for each station, each model, each scenario, over periods 1990-2014 and 2015-2100.
  - CMIP6/hursmax: daily maximum relative humidity CMIP6 data for each station, each model, each scenario, over periods 1985-2014 and 2015-2100.
  - CMIP6/hursmin: daily minimum relative humidity CMIP6 data for each station, each model, each scenario, over periods 1985-2014 and 2015-2100.
  - CMIP6/tasmax: daily maximum temperature CMIP6 data for each station, each model, each scenario, over periods 1985-2014 and 2015-2100.
  - CMIP6/tasmin: daily minimum temperature CMIP6 data for each station, each model, each scenario, over periods 1985-2014 and 2015-2100.
  - CMIP6/LSF: Land-sea fractions for each station and each model, in %.
- OBS: contains observation data, needed to get final figures.
    - OBS/heatwaves: files containing HHW2 thresholds for Gillot, HHW1 and HHW2 frequencies in Gillot observations, for both 1985-2014 and the decade 1990-2000.
    - OBS/HIM: contains daily mean HI for Gillot observations.
    - OBS/HIX: contains daily maximum HI for Gillot observations.
    - OBS/Other files: daily observations during the period 1985-2014 for all stations. Available variables are: tasmax, tasmin, hursmax, hursmin, HIX, HIN, HIM.

## Scripts

Six folders:

- bias_correction: scripts Python to correct biases of CMIP6 tasmin data (example of Southwest Indian domain) with CDF-t. This method can also be applied to the other variables: tasmax, hursmin, hursmax.
    - tasmin_cdft_hist.py: return bias corrected CMIP6 tasmin data in historical period (1985-2014).
    - tasmin_cdft_fut.py: return bias corrected CMIP6 tasmin data in future period (2015-2100).
- calcul_Heat_Index: Heat_Index.py is a script Python to compute NOAA Heat Index from temperature and relative humidity data.
- final_figures:
   -  final_figures/main_paper: this folder contains the four notebooks needed to plot the figures of the main paper.
   -  final_figures/supporting_information: contains the three notebooks needed to plot the figures of the supporting information of the paper.
- hhe: calcul_RL_and_RP_with_GEV.py contains the script to calculate return levels and return periods with a non-stationnary GEV fit.
- hhw_method_1: contains the two scripts to generate HHW1 files.
    - count_hhw1_days.py: compute frequencies of CMIP6 and observation HHW1.
    - find_characteristics_of_hhw1.py: compute characteristics of CMIP6 HHW1 (duration, max intensity).
- hhw_method_2: contains the four scripts to generate HHW2 files.
    - calcul_hhw2_CMIP6.py: calculate properties of HHW2 found in CMIP6 models (max intensity, duration, frequency, severity)
    - calcul_hhw2_OBS.py: calculate properties of HHW2 found in observations (max intensity, duration, frequency, severity)
    - calcul_uncertainties_count_days_hhw2_CMIP6.py: compute confidence interval for CMIP6 HHW2 frequencies.
    - find_characteristics_hhw2_synthesis.py: compute a synthesis of HHW2 characteristics per decade (max intensity, duration, severity).

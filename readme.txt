This is the description of the code and data provided for the JRSS-A manuscript "Pollution State Modeling for Mexico City." All code for our analyses is provided here.  We have not corrected all working directories, so there may be some need to correct these if you use this code. Please feel free to contact Phil White at philawhite@gmail.com if you have questions about the code.


Code/

crps.R  - Defines crps and es. This could be replaced using package ScoringRules
explore.R - Prepares data
pairmean.cpp - Called by crps.R to compute CRPS and ES 

final_model_all_data_emergency_phases/   --- emergency phases scripts
final_mod.R - Fit model to all data, used for emergency phase inference
functions_AR_Final.R - Gibbs sampler and other functions - uses c++ functions 
Gibbs_updates.cpp - Functions for Gibbs sampler
plotting_all.R - script used to create plots
prob_analysis.R - script that takes predictions and infers phase status

final_model_month_exceedance/     --- sequential model fitting
April_fit.R - Fit model through march, predict april - change month number to get different month (e.g. September, fit model through August)
functions_AR_month.R - Gibbs sampler and other functions - uses c++ functions 
Gibbs_updates.cpp - Functions for Gibbs sampler
plotting_april.R - script used to create plots
prob_analysis_april.R - script that takes predictions and infers compliance

heteroscedastic_model_holdout/    ---- Not used (but referred to) in paper
functions_AR_daily_heteroscedastic.R - Sampler that allows variance to change as function of day
Gibbs_updates_het.cpp - Functions for Gibbs sampler
model_heteroscedastic_with_holdout.R - Fits model using heteroscedastic model

homoscedastic_model_holdout/      --- model comparison
functions_AR.R - Sampler for square-root ozone, log-pm10
Gibbs_updates_het.cpp - Functions for Gibbs sampler
model_heteroscedastic_with_holdout.R - Fits model using homoscedastic model


Data/
coord_raw.RData - lat/lon coordinates for stations
O3.csv - 2017 ozone levels
O3_Dec.csv - 2016 December ozone levels
PM10.csv - 2017 pm10 levels
PM10_dec.csv - 2016 December pm10 levels
RH.csv - relative humidity
TMP.csv - temperature in C
UTM_Coordinates - Eastings and Northings for stations

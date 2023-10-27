# PEPRMT-Tidal
PEPRMT-Tidal is a one-dimensional process-based model that predicts gross primary productivity, ecosystem respiration and methane exchange in tidal wetlands at the daily time step.

If you plan to use the model code or data provided in this repository please cite Oikawa et al. In prep.
This repo contains all the files needed to run PEPRMT-Tidal, an updated version of the PEPRMT model.
This model is described in detail in Oikawa et al. In prep.

About the model:
1. Originally PEPRMT (the Peatland Ecosystem Photosynthesis and Methane Transport) model 
was parameterized for restored freshwater wetlands in the Sacramento-San Joaquin River Delta, CA USA
Oikawa et al. 2017 https://doi.org/10.1002/2016JG003438
Presented here is an updated version that works for tidal wetlands (Oikawa et al. In prep)

3. All PEPRMT modules use the same input structure (data) 
however not all models use all variables in the structure;
4. All variables are at the daily time step.
5. Modules are run in succession, first GPP, then Reco and last CH4
6. The GPP module using a light use efficiency equation to predict GPP
    -GPP can be predicted using leaf area index (LAI) or a greenness index from Phenocam data or remote sensing such as EVI or NDVI
    -PEPRMT-Tidal applied in Oikawa et al. In prep uses EVI from Landsat
7. The Reco module uses a Dual Arrhenius Michaelis-Menten (DAMM) approach to predict Reco from 2 carbon pools
8. The CH4 module also uses the DAMM approach and includes inhibition factors to decrease CH4 production in presence of NO3 and SO4

About the data:
We have included here a dataset used in the manuscript which sources data from Ameriflux eddy covariance towers across the United States.
These data have been filtered and gapfilled following methods outlined in Oikawa et al. 2023.
Ancillary data from local tidal streams included continuous NO3 or salinity measurements are also included.


About the Rmd file:
We have also included an Rmd file that loads in the dataset and runs the PEPRMT Tidal modules in sequence.
This Rmd also includes plotting of the time series for each site.
We hope this helps future users to better understand how the data need to be organized and how the modules run in sequence and feed into each other.

If you have any questions, please contact Patty Oikawa patty.oikawa@gmail.com
Thank you!

   

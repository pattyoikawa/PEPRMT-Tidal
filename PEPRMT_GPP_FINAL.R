###################################################
###  PEPRMT-Tidal Gross Primary Productivity (GPP) Module v1.0### 
###################################################

#Developed by Patty Oikawa
#patty.oikawa@gmail.com
#Published in December 2023

#About the model:
#1. Originally PEPRMT (the Peatland Ecosystem Photosynthesis and Methane Transport) model 
  #was parameterized for restored freshwater wetlands in the Sacramento-San Joaquin River Delta, CA USA
  # Oikawa et al. 2017 https://doi.org/10.1002/2016JG003438
  #Presented here is an updated version that works for tidal wetlands (Oikawa et al. 2023) 
#2. All PEPRMT modules use the same input structure (data) 
#   however not all models use all variables in the structure;
#3. All variables are at the daily time step.
#4. Modules are run in succession, first GPP, then Reco and last CH4
#5. This model predicts GPP using a light use efficiency equation
    #GPP can be predicted using leaf area index (LAI) or a greenness index from Phenocam data or remote sensing such as EVI or NDVI
    #PEPRMT-Tidal applied in Oikawa et al. 2023 uses EVI from Landsat

#Inputs: 
#1. Theta: a vector of 4 parameter values that were determined via MCMC Bayesian fitting 
#theta <- c(0.7479271,   1.0497113, 149.4681710,  94.4532674 )
#See Oikawa et al. 2023 

#2. Data: a data frame containing 15 variables as described at start of function.


#Outputs:
# a data frame containing
# 1. GPP:gross primary productivity
# 2. APAR:absorbed photosynthetically active radiation

#units of output
# 1. GPP: g C CO2 m^-2 day^-1
# 2. APAR: #umol m-2 d-1


PEPRMT_GPP_final <-  function(theta,data) {

  
  #---CREATE A SPACE TO COLLECT ITERATIVE RESULTS---#
  q <- unique(as.integer(data$site))
  outcome_lst <- vector('list', length(q))   
  
  #---ITERATIVE LOOP TO RUN THE MODEL ACROSS DIFFERENT SITES---#
  for(i in 1:length(q)){
    #  subset your data here, then create the exogenous variables here
    d <- subset(data, site == i)
    
    #Exogenous Variables
    Time_2 = d[,1] # day of year (1-infinite # of days)
    DOY_disc_2=d[,2] #discontinuous day of year that starts over every year (1-365 or 366)
    Year_2=d[,3] #year 
    TA_2 = d[,4] #Air temperature (C)
    WT_2 = d[,5] #water table depth (cm) equals 0 when water table at soil surface
    PAR_2 <- d[,6] #photosynthetically active radiation (umol m-2 d-1)
    LAI_2 <- d[,7] #Leaf area index (if not using can be 0s or NaN)
    GI_2 <- d[,8] #greeness index from Phenocam (GCC) or Landsat EVI etc (unitless)
    FPAR <- d[,9] #If using LAI data, set FPAR variable to 1's, if using a greenness index set FPAR to 0's
    LUE<- d[,10] #growing season LUE computed for each site using measured GPP in gC per umol
    wetland_age_2= d[,11] #Age of wetland in years
    Sal <- d[,12] #Salinity (ppt)
    NO3 <- d[,13] #Dissolved NO3 (mg/L)
    #Season_drop_2 <- d[,13] #not used in PEPRMT-Tidal (was used in original PEPRMT model in peatlands) 
                      #Season variable that is set to 1 in winter (DOY 1-88, 336-365), 2 pre-spring (DOY 89-175), 3 spring (DOY 176-205), 4 summer (DOY 206-265), 5 fall (DOY 266-335)
    SOM_2 <-d[,14] #Decomposed Organic matter : all the decomposed soil organic matter in top meter of soil informed buy MEM inclusive of current year
    site_2 <-d[,15] #Site: if running more than 1 site, have 1s in this column for first site, 2s for 2nd site and so on
  

##########COMPUTE GPP################################
#PARAMETERS
Ha <- theta[3]+30 #default=30;#activation energy for general crop plant (KJ mol-1)
Hd <- theta[4]+100 #default=100;(KJ mol-1)
#CONSTANTS
vcopt <- 1.0
R_t <- 0.00831 #KJ mol-1 K-1
T_opt <- 25 + 274.15 #(K); our Temp opt for Ps is 25C
  
#EQUATIONS
#Decide how to compute fPAR
#If Fpar = 1 then calculate FPAR using LAI if FPAR =0 then calculate using GI 
if (sum(FPAR)>0) {
  k <- 0.8 #ranges between 0-1
  fPAR_2 <- 0.95*(1-exp(-k*LAI_2)) #for an LAI=4.9, fpar=0.87--Yuan 2007
  fPAR_2 <- (fPAR_2/2)*10^4
} else {
  fPAR_2 <- theta[1]+theta[2]*GI_2
}


APAR_2 <- fPAR_2*PAR_2 #umol m-2 daily average
  
AirT_K <- TA_2 + 274.15 #C to Kelvin


vct <- vector("numeric",length(Time_2))
NPP_FPAR_T <- vector("numeric",length(Time_2))

for (t in 1:length(Time_2)) { 
  exponent1 <- (Ha*(AirT_K[t]-T_opt))/(AirT_K[t]*R_t*T_opt)
  exponent2 <- (Hd*(AirT_K[t]-T_opt))/(AirT_K[t]*R_t*T_opt)
  top <- Hd*exp(exponent1)  
  bottom <- Hd-(Ha*(1-exp(exponent2)))
  vct[t] = vcopt*(top/bottom)

  NPP_FPAR_T[t] <- ((vct[t]*(APAR_2[t]*LUE[t]))) #umol m-2 d-1* gC/umol == g C m-2 d-1


}

GPP <- (NPP_FPAR_T)*-1 #stay as g C m-2 d-1 where negative values= uptake

  w <- cbind(GPP,APAR_2,Time_2) %>%
    as.data.frame(.)
  # store d in a vector  
  outcome_lst[[i]] <- (w) 
  
}

# combine iterations of loop and return all results
GPP_output <- do.call('rbind', outcome_lst) %>% 
  as.data.frame(.)

return(GPP_output)
}



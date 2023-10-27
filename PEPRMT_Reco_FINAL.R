###################################################
###  PEPRMT-Tidal Ecosystem Respiration (Reco) Module v1.0### 
###################################################

#Developed by Patty Oikawa
#patty.oikawa@gmail.com
#Published in November 2023

#About the model:
#1. Originally PEPRMT (the Peatland Ecosystem Photosynthesis and Methane Transport) model 
  #was parameterized for restored freshwater wetlands in the Sacramento-San Joaquin River Delta, CA USA
  # Oikawa et al. 2017 https://doi.org/10.1002/2016JG003438
  #Presented here is an updated version that works for tidal wetlands (Oikawa et al. 2023)
#2. All PEPRMT modules use the same input structure (data) however not all models use all variables in the structure;
#3. All variables are at the daily time step.
#4. Modules are run in succession, first GPP, then Reco and last CH4


#Inputs: 
#1. Theta: a vector of 4 parameter values that were determined via MCMC Bayesian fitting 
#theta <- c(18.41329, 1487.65701,   11.65972,   61.29611 )
#See Oikawa et al. 2023 

#2. Data: a data frame containing 16 variables (columns) as described at start of function.


#3. Wetland_type
# - "wetland_type == 1" corresponds to a "freshwater peatland"
# - "wetland_type == 2" corresponds to a "tidal wetland" 

#Outputs:
# a data frame containing
# 1.Reco_full: total amount of CO2 emitted by ecosystem respiration (including heterotrophic and autotrophic respiration).
# 2.NEE_mod: Net ecosystem exchanged of CO2 (adding Reco with GPP predicted in GPP module)
# 3.S1: amount of carbon stored in soil carbon pool 1, the labile pool
# 4.S2: amount of carbon stored in soil carbon pool 2, the SOC pool

#units of output
# 1.Reco_full: g C CO2 m^-2 day^-1
# 2.NEE_mod: g C CO2 m^-2 day^-1
# 3.S1: g C  m^-3 (includes top m3 of soil)
# 4.S2: g C  m^-3 (includes top m3 of soil)

library(tidyverse)
library(magrittr)

PEPRMT_Reco_FINAL <- function(theta,
                              data,
                              wetland_type) {
  
  #SET UP Reco
  #SOM
  alpha1 <- 3e3 #g C m-2 d-1;--SET AS CONSTANT; 
  ea1 <- theta[1]*1000 #J mol-1
  km1 <- theta[2] # g C m-3

  #LABILE
  alpha2 <- 3e3 #g C m-2 d-1 --SET AS CONSTANT
  ea2 <- theta[3]*1000 #J mol-1
  km2 <- theta[4] # g C m-3
  
  #initialize labile C pool-set to 0 initially
  C2_init <- 0 #in g C m-3
  
  #Empirical parameters used to comput inhibition of Reco when WTD is high
  a1 <- 0.00033
  a2 <- 0.0014
  a3 <- 0.75
  

  
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
    GPP_2 = d[,16] #Modeled GPP - use output from PEPRMT-GPP (gC m-2 day-1) where negative fluxes = uptake
    
    
    # #Static C allocation theme
    NPPsum_avail_2 <- (c(GPP_2)*-1) #g C m-2 day-1 change to + numbers & give Reco access to GPP
    
    #Time Invariant
    R <- 8.314 #J K-1 mol-1
    RT <- R* (TA_2 + 274.15) #T in Kelvin-all units cancel out
    Vmax1 <- alpha1* exp(-ea1/RT) #g C m-2 d-1 SOM
    Vmax2 <- alpha2* exp(-ea2/RT) #g C m-2 d-1 labile
    
    
  #preallocating space
  S1sol <- vector("numeric",length(Time_2))
  S2sol <- vector("numeric",length(Time_2))
  R1 <- vector("numeric",length(Time_2))
  R2 <- vector("numeric",length(Time_2))
  S1 <- vector("numeric",length(Time_2))
  S2 <- vector("numeric",length(Time_2))
  percent_reduction <- vector("numeric",length(Time_2))
  percent_enhancement <- vector("numeric",length(Time_2))
  Reco_1 <- vector("numeric",length(Time_2))
  Reco_full <- vector("numeric",length(Time_2))
  Ps <- vector("numeric",length(Time_2))
  C2in <- vector("numeric",length(Time_2))
  C1_init <- vector("numeric",length(Time_2))
  percent_available <- vector("numeric",length(Time_2))
  
  for (t in 1:length(Time_2)) { 
    #C allocation
    #only 50% of GPP is available 
    C2in[t] <- NPPsum_avail_2[t]*0.5 # gC m-2 d-1

    C1_init[t] <- SOM_2[t]#"decomposed" Organic matter all the soil organic matter in top meter from MEM inclusive of current year
    
    #if (t == 1 | site_change[t]>0 #if beginning of model or switch sites, start C1 pool over
    if (t == 1) {
      S1[t] <- C1_init[t] #substrate avail NOT affected by water avail-- SOM pool
      S2[t] <- C2_init + C2in[t]  # Ps C pool-- some initial Ps C lingering in soil + day 1 GPPavail
    } else {
      S1[t] <- S1sol[t-1]+C1_init[t]
      S2[t] <- S2sol[t-1]+ C2in[t] #substrate availability based on Ps on time step previous
    }
    
    #Empirical factor for increased availability of SOC during the first 3 yrs following restoration
    if (wetland_age_2[t]<1) {
      percent_available[t] <- 0.6
    } else {
      percent_available[t] <- 1 #for peatlands was 20% now 100% is available 
    }
    
    S1[t] = S1[t]*percent_available[t]  #SOM pool
    
    
    #following Davidson and using multiple eq for diff substrate pools
    R1[t] <- Vmax1[t] * S1[t]/(km1 + S1[t]) #g C m2 d-1 Reaction velocity
    R2[t] <- Vmax2[t] * S2[t] /(km2 + S2[t]) #g C m2 d-1   
    if (R1[t]<0) {R1[t]=0}
    
    if (R2[t]<0) {R2[t]=0}
    
    
    #Reco is reduced by 25% when WT is at or above soil surface
    #--McNicol Silver 2015
    #   a1 <- 0.00033
    #   a2 <- 0.0014
    #   a3 <- 0.75
    #   WT_ex <- c(-30, -20, -10, 0)
    #   percent_red_ex <- c(1, 0.85, 0.77, 0.75)
    
    #   plot(percent_red_ex ~ WT_ex)
    
    # insert logic tree here for freshwater peatland or tidal wetland switch
    #Reduce Reco when WTD is high--not applied to Tidal sites
    if(wetland_type == "1"){
      percent_reduction[t] <- (a1*WT_2[t]^2) - (a2*WT_2[t]) + a3
      if (WT_2[t]>5) {percent_reduction[t] <- 0.75}
      if (percent_reduction[t]>1.25) {percent_reduction[t] <- 1.25}
      if (percent_reduction[t]<0.75) {percent_reduction[t] <- 0.75}
      
      R1[t] <- R1[t]*percent_reduction[t] #g C m2 d-1  Reaction velocity
      R2[t] <- R2[t]*percent_reduction[t] #g C m2 d-1 
    } else {}
    
    #Empirical factor for elevated Reco during the first 3 yrs following restoration
    if (wetland_age_2[t]<4) {
      percent_enhancement[t] <- 1.2
    } else {
      percent_enhancement[t] <- 1
    }
    
    R1[t] = R1[t]*percent_enhancement[t] #umol m2 sec Reaction velocity
    R2[t] = R2[t]*percent_enhancement[t] #umol m2 sec
    
    if (t==1) {
      S1sol[t] = C1_init[t] - (R1[t]) #accounts for depletion of C sources in soil due to Reco and methane production
      S2sol[t] = (C2_init+C2in[t]) - (R2[t])
    } else {
      S1sol[t] = S1[t] - R1[t]
      S2sol[t] = S2[t] - R2[t]
    }
    
    if (S1sol[t]<0) {S1sol[t] <- 0}
    if (S2sol[t]<0) {S2sol[t] <- 0}
    
    ###########EDITED OUT IN CURRENT VERSION############
    #in autumn time or season 6, labile PS C pool empties into SOM
    # if (Season_drop_2[t]>5) {
    #   S1sol[t] = S1sol[t]+(0.2*S2sol[t]) #move part of labile C into SOM pool--mimicing plant matter dying
    #   S2sol[t] = S2sol[t]-(0.2*S2sol[t])
    # }
    # 
    # #in winter time or season 1, labile PS C pool empties into SOM
    # if (Season_drop_2[t]<2) {
    #   S1sol[t] = S1sol[t]+(S2sol[t]) #move entire labile C into SOM pool--mimicing plant matter dying
    #   S2sol[t] = 0
    # }
    ########################################
    
    Reco_1[t] <- R1[t] + R2[t] 
    Reco_full[t] <- (R1[t]) + (R2[t]) #umol m2 d-1
  }
  
  NEE_mod <- GPP_2 + Reco_1 #umol m-2 d-1
  
  w <- cbind(Reco_full, NEE_mod, S1, S2) %>%
    as.data.frame(.)
  # store d in a vector  
  outcome_lst[[i]] <- (w) 

  }
  
  # combine iterations of loop and return all results
  Reco_output <- do.call('rbind', outcome_lst) %>% 
    as.data.frame(.)
  
  return(Reco_output)
}
 


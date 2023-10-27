###################################################
###  PEPRMT-Tidal Methane Module v1.0### 
###################################################

#Developed by Patty Oikawa
#patty.oikawa@gmail.com
#Published in November 2023

#About the model:
#1. Originally PEPRMT (the Peatland Ecosystem Photosynthesis and Methane Transport) model 
  #was parameterized for restored freshwater wetlands in the Sacramento-San Joaquin River Delta, CA USA
  # Oikawa et al. 2017 https://doi.org/10.1002/2016JG003438
  #Presented here is an updated version that works for tidal wetlands 
  #and inhibits CH4 production in response to salinity and nitrate (Oikawa et al. 2023)
#2. All PEPRMT modules use the same input structure (data) however not all models use all variables in the structure;
#3. All variables are at the daily time step.
#4. Modules are run in succession, first GPP, then Reco and last CH4



#Inputs: 
#1. Theta: a vector of 7 parameter values that were determined via MCMC Bayesian fitting 
#theta= <- c( 14.9025078, 0.4644174, 16.7845002, 0.4359649, 15.8857612,0.5120464, 486.4106939, 0.1020278 )
#See Oikawa et al. 2023 

#2. Data: a data frame containing 18 variables described at start of function.
 
#3. Wetland_type
# - "wetland_type == 1" corresponds to a "freshwater peatland"
# - "wetland_type == 2" corresponds to a "tidal wetland" 

#Outputs:
# a data frame containing
# 1. pulse_emission_total: total amount of methane emitted, which is the sum of plant-mediated + diffusive water fluxes.
# 2. Plant_flux_net: net amount of methane released from plants.
# 3. Hydro_flux: net amount of methane that transfer from water to atm.
# 4. M1: pool of methane produced from soil carbon pool 1, the labile pool
# 5. M2: pool of methane produced from soil carbon pool 2, the SOC pool
# 6. trans2: fraction of CH4 released via plant-mediated transport

#units of output
# 1. pulse_emission_total: g C methane m^-2 day^-1
# 2. Plant_flux_net: g C methane m^-2 day^-1
# 3. Hydro_flux: g C methane m^-2 day^-1
# 4. M1: g C methane m^-3 (includes top m3 of soil+water)
# 5. M2: g C methane m^-3 (includes top m3 of soil+water)
# 6. trans2: unitless

library(tidyverse)
library(magrittr)

PEPRMT_CH4_FINAL <- function(theta,
                             data,
                             wetland_type){
  #Constants
  R = 8.314 #J K-1 mol-1

  #CH4 PARAMETERS
  #SOC pool
  M_alpha1 = 6.2e13# gC m-3 d-1 
  M_ea1 = (theta[1]+67.1)*1000 #parameter in kJ mol-1 multiplied by 1000 = J mol-1
  M_km1 =theta[2]+17#g C m-3 
  #Labile C pool
  M_alpha2 = 6.2e14# gC m-3 d-1 
  M_ea2 = (theta[3]+71.1)*1000 #J mol-1
  M_km2 =theta[4]+23#g C m-3 
  
  #CH4 oxidation parameters
  M_alpha3 = 6.2e13# gC m-3 d-1 
  M_ea3 = (theta[5]+75.4)*1000 #J mol-1
  M_km3 =theta[6]+23#g C m-3 
  
  #Salinity sulfate parameters
  kiSO4 = theta[7] #mg L^-1 
  kiNO3 = theta[8] #mg L^-1 
  
  #Parameters for hydrodynamic flux
  k=0.04 #gas transfer velocity (m day-1)
  
  #Parameters for plant-mediated transport
  Vtrans=0.24 #gas transfer velocity through plants(m d-1)
  Oxi_factor=0.35 #percent oxidized during transport
  
  #empirical factors for inhibition of ch4 production when WT drops
  beta1=0.48
  beta2=-0.18
  beta3=0.0042
  #empirical factors for decaying inhibition of CH4 production following first
  #flooding of wetland
  zeta1=5.1e-6 #7.4e-6
  zeta2=0.00058
  zeta3=0.11

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
    S1_2 = d[,17] #Modeled SOC pool - use output from PEPRMT Reco model -cumulative grams C m-3
    S2_2 = d[,18] #Modeled labile C pool - use output from PEPRMT Reco model -cumulative grams C m-3
    
    WT_2_adj=(WT_2/100)+1
    GPPmax=max(GPP_2*-1)
    
    #Time Invariant
    RT = R * (TA_2 + 274.15)#T in Kelvin-all units cancel out
    M_Vmax1 = M_alpha1 * exp(-M_ea1/RT)#g C m-2 d-1 
    M_Vmax2 = M_alpha2 * exp(-M_ea2/RT)#gC m-2 d-1 
    M_Vmax3 = M_alpha3 * exp(-M_ea3/RT)#gC m-2 d-1 
    
    #preallocating space
    S1sol = vector('numeric',length(Time_2))
    S2sol = vector('numeric',length(Time_2))
    
    M1 = vector('numeric', length(Time_2))
    M2 = vector('numeric', length(Time_2))
    M1_full = vector('numeric', length(Time_2))
    M2_full = vector('numeric', length(Time_2))
    M_full = vector('numeric', length(Time_2))
    M_percent_reduction= vector('numeric', length(Time_2))
    M_percent_reduction_2= vector('numeric', length(Time_2))
    
    CH4water= vector('numeric', length(Time_2))
    Hydro_flux= vector('numeric', length(Time_2))
    Plant_flux= vector('numeric', length(Time_2))
    Plant_flux_net= vector('numeric', length(Time_2))
    CH4water_store= vector('numeric', length(Time_2))
    CH4water_0= vector('numeric', length(Time_2))
    Oxi_full= vector('numeric', length(Time_2))
    R_Oxi= vector('numeric', length(Time_2))
    CH4water_0_2= vector('numeric', length(Time_2))
    #Vtrans=zeros(1,length(Time_2))
    #Oxi_factor=zeros(1,length(Time_2))
    trans2=vector('numeric', length(Time_2))
    S1= vector('numeric', length(Time_2))
    S2= vector('numeric', length(Time_2))
    conc_so4AV= vector('numeric', length(Time_2))
    
  #--METHANE TRANSPORT ACROSS DATA---
  for(t in 1:length(Time_2)) {
    
    #parameter for plant-mediated transport--function of GPP
    trans2[t]=((GPP_2[t]+(GPPmax))/GPPmax)
    if (trans2[t]<0) {trans2[t]=0}
    if (trans2[t]>1) {trans2[t]=1}
    
    conc_so4AV[t]=7.4e-8*(Sal[t]*1e6)*1000 # ppm or mg L-1

    M1[t] = (M_Vmax1[t] * (S1_2[t] /(M_km1 + S1_2[t]))) * (kiSO4 /(kiSO4 + conc_so4AV[t])) * (kiNO3 /(kiNO3 + NO3[t]))  #gC m2 d-1 
    
    M2[t] = (M_Vmax2[t] * (S2_2[t] /(M_km2 + S2_2[t])))  * (kiSO4 /(kiSO4 + conc_so4AV[t]))* (kiNO3 /(kiNO3 + NO3[t]))   #gC m2 d-1
    
    if (M1[t]<0) { M1[t]=0}
    if (M2[t]<0) { M2[t]=0}
    
    # For Freshwater peatlands ONLY
    #Empirical eq for CH4 inhibition when WT falls below soil surface
    # if WT below soil surface any time in 10days previous, CH4 production is reduced
    if(wetland_type == 1){
      if (t<=20 & WT_2_adj[t]<0.9){ 
        M_percent_reduction[t]=(beta1*WT_2_adj[t]^2)+(beta2*WT_2_adj[t])+beta3
      } else {M_percent_reduction[t]=1}
      
      if (t>20){ Sel=WT_2_adj[(t-19):t] } else {Sel=5}
      
      if (t>20 & min(Sel, na.rm=T)<0.9) { 
        M_percent_reduction[t]=(beta1*WT_2_adj[t]^2)+(beta2*WT_2_adj[t])+beta3
      } else if (t>20){ M_percent_reduction[t]=1}
      
      if (WT_2_adj[t]<0){  M_percent_reduction[t]=0}
      
      if (M_percent_reduction[t]<0){  M_percent_reduction[t]=0}
      if (M_percent_reduction[t]>1){  M_percent_reduction[t]=1}
      
      #Empirical eq Oikawa for CH4 inhibition for 1 yr following restoration
      if (wetland_age_2[t]<2){
        M_percent_reduction_2[t]=(zeta1*DOY_disc_2[t]^2)+(zeta2*DOY_disc_2[t])+zeta3
      } else { M_percent_reduction_2[t]=1 }
      
      if (M_percent_reduction_2[t]>1) { M_percent_reduction_2[t]=1 }
      if (M_percent_reduction[t] < 0.75){M_percent_reduction[t] = 0.75}
      
      M1[t] = M1[t]*M_percent_reduction[t]  #gC m2 d Reaction velocity
      M2[t] = M2[t]*M_percent_reduction[t]  #gC m2 d
      M1[t] = M1[t]*M_percent_reduction_2[t] #gC m2 d Reaction velocity
      M2[t] = M2[t]*M_percent_reduction_2[t] #gC m2 d
    } else {}
    
    #S1sol and S2sol are the new SOC and labile pools adjusted for C lost through CH4 production
    S1sol[t] = S1[t] - (M1[t]) 
    S2sol[t] = S2[t] - (M2[t])
    
    if (S1sol[t]<0) { S1sol[t]=0}
    if (S2sol[t]<0) { S2sol[t]=0}
    
    M_full[t]=(M1[t]+M2[t])#*fSal #total CH4 produced at this time step in gC m-3 soil day-1
    
    #---COMPUTE CH4 TRANSPORT---
    #make sure WT_2_adj is never negative
    if (WT_2_adj[t]<0) { WT_2_adj[t]=0 }
    
    #now start ch4 transport loop
    if (t==1){
      
      # WT_2_adj = 1 at soil surface
      if (WT_2_adj[t]>1) {
        
        #Methane oxidation is zero when WT above surface
        R_Oxi[t] = 0
        Oxi_full[t]=0
        
        #This assumes you start out the year with no CH4 in water
        #where CH4water_0 is the initial concentration of CH4 in water
        CH4water_0[t]=0 # 
        CH4water_0_2[t]=0 # CH4 produced in previous time step
        #Only modeling CH4 dissolved in the water that goes down 1 m3 into soil
        CH4water[t]= ((M_full[t]*1)+(CH4water_0[t]*WT_2_adj[t]))/ WT_2_adj[t]#gC per m^3 - concentration in water doesn't change
 
        #based on the concentrations in the soil and water, calculate hydro and plant-mediated fluxes
        Hydro_flux[t]=k*CH4water[t]  #gC CH4 m-2 day-1  Hydrodynamic flux 
        Plant_flux[t]=(Vtrans*CH4water[t])*trans2[t]  ##gC CH4 m-2 day-1  Bulk Plant mediated transport 
        Plant_flux_net[t]=Plant_flux[t]*Oxi_factor  ##gC CH4 m-2 day-1  Plant mediated transport after oxidation
        
        #subtract the moles of methane lost from the pools (soil and water) to the atm
        CH4water_store[t]=CH4water[t]-Hydro_flux[t]-Plant_flux[t]  #gC CH4 m-3 stored in the system
        
        #If you start out the year with no water above the surface      
      } else {
        
        # Gives CH4 concentration in water which is now less that 1 m3
        CH4water_0[t]=(M_full[t]*1)+(0.00001*WT_2_adj[t])/WT_2_adj[t]#gC CH4 m-3 - concentration in soil and water are the same
        
        # Methane oxidation turns on when WT falls below the surface
        R_Oxi[t] = M_Vmax3[t] * CH4water_0[t] /(M_km3 + CH4water_0[t])  #gC CH4 m-3 Reaction velocity
        Oxi_full[t]=R_Oxi[t] # gC m-3 d-1
        CH4water[t]=CH4water_0[t]-Oxi_full[t] 
        if (CH4water[t]<0) {  CH4water[t]=0 }
        
        CH4water_0_2[t]=0
        
        #this hydroflux uses a 2nd k parameter which basically inhibits
        #diffusive flux when there is no water above the soil surface
        Hydro_flux[t]=k*CH4water[t]  #gC m-3 d-1 Hydrodynamic flux 
        Plant_flux[t]=(Vtrans*CH4water[t])*trans2[t]  #gC m-3 d-1 Bulk Plant mediated transport 
        Plant_flux_net[t]=Plant_flux[t]*Oxi_factor  #gC m-3 d-1 Plant mediated transport after oxidation
        
        CH4water_store[t]=CH4water[t]-Plant_flux[t]-Hydro_flux[t]  #gC m-3 stored in the system
      }
      
      # If you have water, then the CH4 should mix between the 2 layers and concentrations should be the same in water and soil
    } else {
      
      if (WT_2_adj[t]>1) {
        
        #account for any changes in concentration of CH4 due to any change in WT_2_adj height
        CH4water_0[t]=(CH4water_store[t-1] * WT_2_adj[t-1]) / WT_2_adj[t]
        
        #Methane oxidation is zero when WT above surface
        R_Oxi[t] = 0
        Oxi_full[t]=0 
        CH4water_0_2[t]=0
        #Now add the new CH4 produced today to the soil and let it increase in concentration
        CH4water[t]= ((M_full[t]*1)+(CH4water_0[t]*WT_2_adj[t]))/ WT_2_adj[t]#gC m-3 - concentration in water doesn't change
        
        #again compute fluxes
        Hydro_flux[t]=k*CH4water[t]  #gC m-3 d-1 Hydrodynamic flux 
        Plant_flux[t]=(Vtrans*CH4water[t])*trans2[t]  #gC m-3 d-1 Plant mediated transport 
        Plant_flux_net[t]=Plant_flux[t]*Oxi_factor #gC m-3 d-1 Plant mediated transport after oxidation
        
        #subtract the moles of methane lost from the pools (soil and #water) to the atm
        CH4water_store[t]=CH4water[t]-Plant_flux[t]-Hydro_flux[t] #gC m-3 stored in the system
        
      } else {
        #if you don't have WT_2_adj above the soil, then all the CH4 in the water goes to the soil
        #First, account for increased concentration of CH4 in soil now that WT_2_adj has dropped
        CH4water_0[t]=(CH4water_store[t-1] * WT_2_adj[t-1]) / WT_2_adj[t]
        
        #now add new CH4 to this new concentrated pool of CH4 in the soil
        CH4water_0_2[t]= ((M_full[t]*1)+(CH4water_0[t]*WT_2_adj[t]))/ WT_2_adj[t]#gc per m^3 - concentration in water doesn't change
        
        #Methane oxidation turns on when WT falls below the surface
        R_Oxi[t] = M_Vmax3[t] * CH4water_0_2[t] /(M_km3 + CH4water_0_2[t])  #gC m2 d Reaction velocity
        Oxi_full[t]=R_Oxi[t]#gC m-3 day-1
        CH4water[t]=CH4water_0_2[t]-Oxi_full[t]#now you have less ch4
        if (CH4water[t]<0) {CH4water[t]=0}
        
        
        Hydro_flux[t]=k*CH4water[t]  #gC m-3 d-1 Hydrodynamic flux 
        Plant_flux[t]=(Vtrans*CH4water[t])*trans2[t]  #gC m-3 d-1 Plant mediated transport
        Plant_flux_net[t]=Plant_flux[t]*Oxi_factor  #gC m-3 d-1 Plant mediated transport after oxidation
        
        CH4water_store[t]=CH4water[t]-Plant_flux[t]-Hydro_flux[t] #gC m-3 stored in the system
      }
    }
  }
  
  pulse_emission_total = Plant_flux_net+Hydro_flux  #gC CH4 m-2 day-1 total CH4 flux to atm
  
  w <- (data.frame(cbind(pulse_emission_total,
                           Plant_flux_net,
                           Hydro_flux,
                           M1,
                           M2,
                           trans2)))
  
  # store d in a vector  
  outcome_lst[[i]] <- (w) 
  # return(ch4_output)
  }
  
  # combine iterations of loop and return all results
  peprmtCH4 <- do.call('rbind', outcome_lst) %>% 
    as.data.frame(.)
  
  return(peprmtCH4)
}
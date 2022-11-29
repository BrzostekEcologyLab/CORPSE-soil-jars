
#Baseline Input Parameters
## Data frame with definitions of expected parameters 
params_definitions <- data.frame( 
  "Vmaxref" = as.character ('Relative maximum enzymatic decom rates (Fast, Slow, Necro)'), 
  "Ea" = as.character ('Activation engery (Fast, Slow, Necro)'),
  "kC" = as.character ('Michealis-Menton parameter (Fast, Slow, Necro)'),
  "gas_diffusion_exp" = as.character ('Determines suppression of decomp at high soil moisture'), 
  "minMicrobeC" = as.character ('Minimum microbial biomass (fraction of total C)'), 
  "Tmic"= as.character ('Microbial lifetime at 20C (years)'),
  "et" = as.character ('Fraction of turnover not converted to CO2'),
  "eup" = as.character ('Carbon uptake efficiency (Fast, Slow, Necro)'),
  "nup" = as.character ('Nitrogen uptake efficiency (Fast, Slow, Necro)'),
  "tProtected" = as.character ('Protected C turnover time (years)'), 
  "frac_N_turnover_min" = as.character ('Fraction of microbial biomass N turnover that is mineralized'), 
  "protection_rate" = as.character ('Protected carbon formation rate (year-1)'), 
  "CN_Microbe" = as.character ('C:N ratio of microbial biomass for AM sites'), 
  "max_immobilization_rate" = as.character ('Maximum N immobilization rate (fraction per day)'),
  "substrate_diffusion_exp" = as.character ('Determines suppression of decomp at low soil moisture'),
  "frac_turnover_slow" = as.character('Fraction of microbial biomass N turnover that goes to slow pool'),
  "new_resp_units" = as.character ('If TRUE, Vmaxref has units of 1/years and assumes optimal soil moisture has a relative rate of 1.0'),
  "iN_loss_rate" = as.character ('Loss rate of inorganic N pool (year-1) > 1 because it takes much less than a year for it to be removed'))

##Parameter sets: 
##OG: Baseline values
params <- data.frame( "Vmaxref_Fast"= 33, "Vmaxref_Slow"= 0.35, "Vmaxref_Necro"= 22, "Ea_Fast"= 5e3 , "Ea_Slow"= 25e3, "Ea_Necro"= 3e3, "kC_Fast"= 0.007, "kC_Slow"= 0.005, "kC_Necro"= 0.007, "gas_diffusion_exp"= 0.6, "minMicrobeC"= 1e-3, 
                      "Tmic"= 0.25, "et"= 0.6, 
                      "eup_Fast"= 0.6, "eup_Slow"= 0.001, 
                      "eup_Necro"= 0.6, "tProtected"= 100.0,  "frac_N_turnover_min"= 0.2,
                      "protection_rate_Fast"= 0.6, 
                      "protection_rate_Slow"= 0.001,  
                      "protection_rate_Necro"= 4, "nup_Fast"= 0.3, "nup_Slow"= 0.3,  "nup_Necro"= 0.3, "CN_Microbe"= 6.1, "max_immobilization_rate"= 3.65,  "substrate_diffusion_exp"= 1.5, "new_resp_units"= TRUE, "iN_loss_rate"= 5.0,  "frac_turnover_slow"= 0.2)

##Data based:
##Faster turnover (tmic = 1 week or 7/365), protection rate increases x 4, CUE Scaled by data to average = 0.6
paramsNEW <- data.frame( "Vmaxref_Fast"= 33, "Vmaxref_Slow"= 0.35, "Vmaxref_Necro"= 22, "Ea_Fast"= 5e3 , "Ea_Slow"= 25e3, "Ea_Necro"= 3e3, "kC_Fast"= 0.007, "kC_Slow"= 0.005, "kC_Necro"= 0.007, "gas_diffusion_exp"= 0.6, "minMicrobeC"= 1e-3,
                      "Tmic"= 7/365, "et"= 0.6,
                      "eup_Fast"= c(.396,.305,.333,.171)*.6/(sum(.396,.305,.333,.171)/4), "eup_Slow"= 0.001,
                      "eup_Necro"= c(.396,.305,.333,.171)*.6/(sum(.396,.305,.333,.171)/4), "tProtected"= 100.0,  "frac_N_turnover_min"= 0.2,
                      "protection_rate_Fast"= 2.4,
                      "protection_rate_Slow"= 0.004,
                      "protection_rate_Necro"=16, "nup_Fast"= 0.3, "nup_Slow"= 0.3,  "nup_Necro"= 0.3, "CN_Microbe"= 6.1, "max_immobilization_rate"= 3.65,  "substrate_diffusion_exp"= 1.5, "new_resp_units"= TRUE, "iN_loss_rate"= 5.0,  "frac_turnover_slow"= 0.2)

if (NewParams==TRUE){ params=paramsNEW } 

params_bulk <- data.frame( "Vmaxref_Fast"=  33, "Vmaxref_Slow"=  0.6, "Vmaxref_Necro"= 22, "Ea_Fast"= 5e3 , "Ea_Slow"=30e3, "Ea_Necro"= 3e3, "kC_Fast"= 0.01, "kC_Slow"= 0.01, "kC_Necro"= 0.01, "gas_diffusion_exp"= 0.6, "minMicrobeC"= 1e-3, "Tmic"= 0.25, "et"= 0.6, "eup_Fast"= 0.6 , "eup_Slow"=0.1 , "eup_Necro"= 0.6, "tProtected"= 100.0,  "frac_N_turnover_min"= 0.2,
                           "protection_rate_Fast"= 0.6, "protection_rate_Slow"= 0.001,  "protection_rate_Necro"= 4, "nup_Fast"= 0.3, "nup_Slow"= 0.3,  "nup_Necro"= 0.3, "CN_Microbe"= 6.1, "max_immobilization_rate"= 3.65,  "substrate_diffusion_exp"= 1.5, "new_resp_units"= TRUE, "iN_loss_rate"= 5.0,  "frac_turnover_slow"= 0.2)

## CORPSE leaf and root litter parameters
chem_types<-c('Fast','Slow','Necro')
claymod <- 1.0  #Scalar that modifies the ability of clays to sorb and protect C 
inorg_Ndep <- 0.01 # inorganic N deposition. Without this the microbes crash and C accumulates infinitely

##Protection rate parameter based on percent clay content from Mayes et al (2012) Table 1
prot_clay <- function (claypercent, slope=0.4833, intercept=2.3282, BD=1.15, porosity=0.4) {
  prot<-1.0*(10**(slope*log10(claypercent)+intercept)*BD*1e-6)
  return (prot)
}

##Calculate rates of change for ALL CORPSE pools 
##T = Temperature (Kelvin)
##Theta = soil water content (fraction of saturation)

##Next 3 functions are needed to run model code
##Functions are ordered to allow code to run

##Function to calculate Vmax of microbial decomposition 
##Vmax function, normalized to Tref=293.15 (T is in Kelvin)
Vmax <- function (T, params, Tref = 293.15, Rugas = 8.314472) {
  Fast <- params$Vmaxref_Fast*exp(-params$Ea_Fast*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))
  Slow <- params$Vmaxref_Slow*exp(-params$Ea_Slow*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))
  Necro <- params$Vmaxref_Necro*exp(-params$Ea_Necro*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))
  Vmax<-data.frame(Fast,Slow,Necro)
  return(Vmax)
}

##Function to sum carbon types  
sumCtypes <- function(SOM, prefix, suffix='C') {
  out <- SOM[[paste(prefix, chem_types[1], suffix, sep='')]]
  if (length(chem_types)>1){
    for (i in 2:length(chem_types)) {
      out <- out+SOM[[paste(prefix,chem_types[i],suffix,sep='')]]
      
    }
  }
  return(out)
}


##Function to calculate Decomposition rate 
decompRate <- function(SOM, T, theta, params){
  theta[theta<0] <- 0.0
  theta[theta>1.0] <- 1.0
  if (params$new_resp_units[1]==TRUE) {
    theta_resp_max<-params$substrate_diffusion_exp[1]/(params$gas_diffusion_exp[1]*(1.0+params$substrate_diffusion_exp[1]/params$gas_diffusion_exp[1]))
    aerobic_max<-theta_resp_max**params$substrate_diffusion_exp[1]*(1.0-theta_resp_max)**params$gas_diffusion_exp[1]
  } else {
    aerobic_max<-1.0
  }
  eq_1 <- theta**params$substrate_diffusion_exp
  eq_2 <- (1.0-theta)**params$gas_diffusion_exp/aerobic_max
  
  vmax <- Vmax(T,params)
  decompRate <- data.frame(matrix(nrow=nspp,ncol=0))
  for (ctypes in chem_types) {
    #drate <- vmax[[ctypes]]*eq_1*eq_2*(SOM[[paste('u',ctypes,'C',sep='')]]*SOM$livingMicrobeC/((sumCtypes(SOM,'u'))*params[[paste('kC',ctypes,sep='_')]]+SOM$livingMicrobeC))
    drate <- vmax[[ctypes]]*(SOM[[paste('u',ctypes,'C',sep='')]]*SOM$livingMicrobeC/((sumCtypes(SOM,'u'))*params[[paste('kC',ctypes,sep='_')]]+SOM$livingMicrobeC))
    zeroMB = SOM$livingMicrobeC==0.0
    drate[zeroMB] <- 0.0
    decompRate[paste(ctypes,'C',sep='')] <- drate
    NC_ratio <- SOM[[paste('u',ctypes,'N',sep='')]]/SOM[[paste('u',ctypes,'C', sep='')]]
    NC_ratio[SOM[[paste('u',ctypes,'C', sep='')]]==0.0]<-0.0
    decompRate[paste(ctypes,'N',sep='')] <- drate*NC_ratio
  }
  return (decompRate)
}


CORPSE_decomp <- function(SOM,T,theta,params,claymod=1.0) {
  ##Calculate maximum potential C decomposition rate
  decomp<-decompRate(SOM,T,theta,params)

  ##Microbial Turnover
  microbeTurnover<-(SOM$livingMicrobeC-params$minMicrobeC*sumCtypes(SOM,'u','C'))/params$Tmic ##kg/m2/yr
  microbeTurnover[microbeTurnover<0.0]<-0.0

  ##Calculating maintenance respiration and creating zeros matrix for overflow_resp
  maintenance_resp <- microbeTurnover*(1.0-params$et)
  overflow_resp <- maintenance_resp*0

  ##Calculating fraction dead microbes in C and N
  deadmic_C_production <- microbeTurnover*params$et ##actual fraction of microbial turnover
  deadmic_N_production<-(microbeTurnover*params$et)/(params$CN_Microbe)

  ##C and N available for microbial growth
  carbon_supply<-vector(mode='numeric', length=nspp)
  nitrogen_supply<-vector(mode='numeric', length=nspp)
  for (ctypes in chem_types) {
    carbon_supply<-carbon_supply+decomp[[paste(ctypes,'C',sep='')]]*params[[paste('eup',ctypes,sep='_')]]
    nitrogen_supply<-nitrogen_supply+decomp[[paste(ctypes,'N',sep='')]]*params[[paste('nup',ctypes,sep='_')]]
  }

  IMM_N_max <- params$max_immobilization_rate*365*SOM$inorganicN/(SOM$inorganicN+params$max_immobilization_rate)
  
  dmicrobeC <- vector(mode='numeric', length=nspp)
  dmicrobeN <- vector(mode='numeric', length=nspp)
  CN_imbalance_term <- vector(mode='numeric', length=nspp)

  ##Growth is nitrogen limited, with not enough mineral N to support it with max immobilization
  ##Model isn't able to do multiple true/false runs at the same time, so localized N limitation, N immonbilization, or C limitation is
  ##individually calculated for each litter type. 
  
  for (j in 1:length(nspp)){##runs each litter type 
    ##1. Determine if loc_Nlim is true: 
  loc_Nlim <- (carbon_supply[j] - maintenance_resp[j] )>((nitrogen_supply[j]  + IMM_N_max[j] ) * params$CN_Microbe[j] )
    ##If true, record data
  if(loc_Nlim==TRUE){
  CN_imbalance_term[j] <- (-IMM_N_max[j])

  dmicrobeC[j] <- ((nitrogen_supply[j] + IMM_N_max[j]) * params$CN_Microbe[j] - microbeTurnover[j] * params$et[j])

  dmicrobeN[j] <- (nitrogen_supply[j] + IMM_N_max[j] - microbeTurnover[j] * params$et[j]/params$CN_Microbe[j])

  overflow_resp[j] <- carbon_supply[j] - maintenance_resp[j] - (nitrogen_supply[j] + IMM_N_max[j]) * params$CN_Microbe[j]
}
  ##Growth must be supported by immobilization of some mineral nitrogen, but is ultimately carbon limited
  ##2. Determine if loc_immob is true
  loc_immob <- (carbon_supply[j] - maintenance_resp[j] >= nitrogen_supply[j]*params$CN_Microbe[j]) & (carbon_supply[j] - maintenance_resp[j] < (nitrogen_supply[j]+IMM_N_max[j])*params$CN_Microbe[j])
  ##if true, record data
  if(loc_immob==TRUE){
  CN_imbalance_term[j] <- (-((carbon_supply[j] - maintenance_resp[j])/params$CN_Microbe[j] - nitrogen_supply[j]))
  
  dmicrobeC[j] <- (carbon_supply[j] - microbeTurnover[j])

  dmicrobeN[j] <- ((carbon_supply[j]-maintenance_resp[j])/params$CN_Microbe[j] - microbeTurnover[j]*params$et[j]/params$CN_Microbe[j])
}
  ##Growth is carbon limited and extra N is mineralized
  ##determine if loc_Clim is true
  loc_Clim <- !(loc_Nlim | loc_immob)
  ##if true,record data
 if(loc_Clim==TRUE){
  dmicrobeC[j] <- (carbon_supply[j] - microbeTurnover[j])

  dmicrobeN[j] <- ((carbon_supply[j] - maintenance_resp[j])/params$CN_Microbe[j] - microbeTurnover[j]*params$et[j]/params$CN_Microbe[j])

  CN_imbalance_term[j] <- nitrogen_supply[j] - (carbon_supply[j]-maintenance_resp[j])/params$CN_Microbe[j]
}
  }
  ##CO2 production and cumulative CO2 produced by cohort
  CO2prod <- maintenance_resp + overflow_resp
  for (ctypes in chem_types) {
    CO2prod <- CO2prod + decomp[[paste(ctypes,'C',sep='')]]*(1.0-params[[paste('eup',ctypes,sep='_')]])
  }

  ##Update protected carbon
  protectedturnover <- data.frame(matrix(nrow=nspp,ncol=0))
  protectedprod <- data.frame(matrix(nrow=nspp,ncol=0))

  for (ctypes in chem_types) {
    protectedturnover[paste(ctypes,'C',sep='')] <- SOM[paste('p',ctypes,'C',sep='')]/params$tProtected

    protectedprod[paste(ctypes,'C',sep='')] <- SOM[paste('u',ctypes,'C',sep='')]*params[[paste('protection_rate',ctypes,sep='_')]]*claymod

    protectedturnover[paste(ctypes,'N',sep='')] <- SOM[paste('p',ctypes,'N',sep='')]/params$tProtected

    protectedprod[paste(ctypes,'N',sep='')] <- SOM[paste('u',ctypes,'N',sep='')]*params[[paste('protection_rate',ctypes,sep='_')]]*claymod
  }

  derivs <- SOM

  derivs$livingMicrobeC <- dmicrobeC
  derivs$livingMicrobeN <- dmicrobeN
  derivs$CO2 <- CO2prod
  derivs$inorganicN <- CN_imbalance_term

  for (ctypes in chem_types) {
    derivs[[paste('u',ctypes,'C',sep='')]] <- (-decomp[[paste(ctypes,'C', sep='')]]) + protectedturnover[[paste(ctypes,'C', sep='')]] - protectedprod[[paste(ctypes,'C', sep='')]]

    derivs[[paste('p',ctypes,'C',sep='')]] <- protectedprod[[paste(ctypes,'C', sep='')]] - protectedturnover[[paste(ctypes,'C', sep='')]]

    derivs[[paste('u',ctypes,'N',sep='')]] <- (-decomp[[paste(ctypes,'N',sep='')]]) + protectedturnover[[paste(ctypes,'N', sep='')]] - protectedprod[[paste(ctypes,'N', sep='')]]

    derivs[[paste('p',ctypes,'N',sep='')]] <- protectedprod[[paste(ctypes,'N', sep='')]] - protectedturnover[[paste(ctypes,'N', sep='')]]
  }

  derivs['uNecroC'] <- derivs['uNecroC'] + deadmic_C_production * (1.0-params$frac_turnover_slow)

  derivs['uSlowC'] <- derivs['uSlowC'] + deadmic_C_production * params$frac_turnover_slow

  turnover_N_min <- deadmic_N_production * params$frac_N_turnover_min

  turnover_N_slow <- deadmic_N_production * params$frac_turnover_slow

  derivs['uNecroN'] <- derivs['uNecroN'] + deadmic_N_production - turnover_N_min-turnover_N_slow

  derivs['uSlowN'] <- derivs['uSlowN'] + turnover_N_slow

  derivs['inorganicN'] <- derivs['inorganicN']+turnover_N_min

  return(derivs)
}



##{Main Code Function: This is what has been edited to allow multiple model runs}

##Function Definition: 
##Inputs have been separated out 
  CORPSE_MAIN <- function(nyears,nspp,litter_added,CORPSE_full_spinup_bulk,soil_T,soil_VWC){
    
    
  #Function Variables:
  timestep <- nyears*365
  CORPSEstep <- 1/365 ## this is daily values
  claymod <- 1.0  #Scalar that modifies the ability of clays to sorb and protect C 
  inorg_Ndep <- 0.01 # inorganic N deposition. Without this the microbes crash and C accumulates infinitely
  
  #Input files: read in from excel, these are defined as model inputs
  # adding litter C in litterbag
  litterbag_int <-read.csv(litter_added)
  #Initial Data
  bulk_int <- read.csv(CORPSE_full_spinup_bulk)
  
  ##soil abiotic conditions (temperature, water content)
  soilT <- as.matrix(read.csv(soil_T))
  soilVWC <- as.matrix(read.csv(soil_VWC))

  
  ## Make an empty data frame with column names for pools (unprotected and protected) and chem_types (Fast, Slow, Necro) for each soil compartment ####WHAT ARE THESE UNITS? ERB
  
  litterbag <- 
    data.frame("uFastC"= numeric (nspp),"uSlowC"= numeric(nspp), "uNecroC"= numeric(nspp), "pFastC"= numeric(nspp), "pSlowC"= numeric(nspp),"pNecroC"= numeric(nspp), "livingMicrobeC"= numeric(nspp),"uFastN"= numeric(nspp), "uSlowN"= numeric(nspp), "uNecroN"= numeric(nspp), "pFastN"= numeric(nspp),"pSlowN"= numeric(nspp), "pNecroN"= numeric(nspp),"inorganicN"= numeric(nspp), "CO2"= numeric(nspp), "livingMicrobeN"= numeric(nspp)
    )
  
  bulk <- 
    data.frame("uFastC"= numeric (nspp),"uSlowC"= numeric(nspp), "uNecroC"= numeric(nspp), "pFastC"= numeric(nspp), "pSlowC"= numeric(nspp),"pNecroC"= numeric(nspp), "livingMicrobeC"= numeric(nspp),"uFastN"= numeric(nspp), "uSlowN"= numeric(nspp), "uNecroN"= numeric(nspp), "pFastN"= numeric(nspp),"pSlowN"= numeric(nspp), "pNecroN"= numeric(nspp),"inorganicN"= numeric(nspp), "CO2"= numeric(nspp), "livingMicrobeN"= numeric(nspp)
    )
  
  # load initial data into dataframes made in step above
  litterbag[1:nspp,1:16] <- litterbag_int[1:nspp,1:16] 
  litterbag$livingMicrobeN <- litterbag$livingMicrobeC/params$CN_Microbe
  
  bulk[1:nspp,1:16] <- bulk_int[1:nspp,1:16] 
  bulk$livingMicrobeN <- bulk$livingMicrobeC/params_bulk$CN_Microbe
  

  
  ##Set up empty matrices to hold CORPSE model outputs 
  
  litterbag_uFastC <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_uSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_uNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_pFastC <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_pSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_pNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_uFastN <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_uSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_uNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_pFastN <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_pSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_pNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_livingMicrobeC <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_inorganicN <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_CO2 <- matrix(NA, nrow=timestep, ncol=nspp)
  litterbag_livingMicrobeN <- matrix(NA,nrow=timestep, ncol=nspp)
  
  bulk_uFastC <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_uSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_uNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_pFastC <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_pSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_pNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_uFastN <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_uSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_uNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_pFastN <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_pSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_pNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_livingMicrobeC <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_inorganicN <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_CO2 <- matrix(NA, nrow=timestep, ncol=nspp)
  bulk_livingMicrobeN <- matrix(NA,nrow=timestep, ncol=nspp)
  
 
  # initializing shared inorganic N pool
  shared_inorganicN <- numeric(nspp)
  
  ## Part 4 Call CORPSE code
  
  for (i in 1:timestep) {
    ##Start CORPSE model main loop
    print(i)
    #Parse step into DOY 
    ##this takes the timestep and makes it into DOY
    k<-(i%%365)
    if(k==0){print(k)}
    if(k==0){k<-365}
    
    ##Get Temperature and theta (soil moisture) values for this time point
    T_step <- soilT[k,1]+273.15
    #porosity <- 0.5
    theta_step <- soilVWC[k,1]#porosity 
    
    ##Running the CORPSE function for LIDET data
  
    #Remove n lim for incubations
    shared_inorganicN<- .2
    
    # Litterbag layer
    litterbag$inorganicN <- shared_inorganicN
    
    results_litterbag <- CORPSE_decomp(litterbag, T_step, theta_step, params, claymod) # units of kg change per year
    

    # Bulk
    bulk$inorganicN <- shared_inorganicN
    
    results_bulk <- CORPSE_decomp(bulk, T_step, theta_step, params_bulk, claymod) # units of kg change per year

    ## Update the pools in SOM by add derivs*dt (length of time step) to each SOM pool.  
    ## This simply converts the units from mass per year to mass per the selected time step
    litterbag <- litterbag + results_litterbag*CORPSEstep  
    litterbag$inorganicN <- shared_inorganicN
    
    bulk <- bulk + results_bulk*CORPSEstep 
    bulk$inorganicN <- shared_inorganicN

    ##Save all pool values in output matrices
    litterbag_uFastC[i,]<-litterbag$uFastC
    litterbag_uSlowC[i,]<-litterbag$uSlowC
    litterbag_uNecroC[i,]<-litterbag$uNecroC
    litterbag_pFastC[i,]<-litterbag$pFastC
    litterbag_pSlowC[i,]<-litterbag$pSlowC
    litterbag_pNecroC[i,]<-litterbag$pNecroC
    litterbag_uFastN[i,]<-litterbag$uFastN
    litterbag_uSlowN[i,]<-litterbag$uSlowN
    litterbag_uNecroN[i,]<-litterbag$uNecroN
    litterbag_pFastN[i,]<-litterbag$pFastN
    litterbag_pSlowN[i,]<-litterbag$pSlowN
    litterbag_pNecroN[i,]<-litterbag$pNecroN
    litterbag_livingMicrobeC[i,]<-litterbag$livingMicrobeC
    litterbag_inorganicN[i,]<-litterbag$inorganicN
    litterbag_CO2[i,]<-litterbag$CO2
    litterbag_livingMicrobeN[i,]<-litterbag$livingMicrobeN
    bulk_uFastC[i,]<-bulk$uFastC
    bulk_uSlowC[i,]<-bulk$uSlowC
    bulk_uNecroC[i,]<-bulk$uNecroC
    bulk_pFastC[i,]<-bulk$pFastC
    bulk_pSlowC[i,]<-bulk$pSlowC
    bulk_pNecroC[i,]<-bulk$pNecroC
    bulk_uFastN[i,]<-bulk$uFastN
    bulk_uSlowN[i,]<-bulk$uSlowN
    bulk_uNecroN[i,]<-bulk$uNecroN
    bulk_pFastN[i,]<-bulk$pFastN
    bulk_pSlowN[i,]<-bulk$pSlowN
    bulk_pNecroN[i,]<-bulk$pNecroN
    bulk_livingMicrobeC[i,]<-bulk$livingMicrobeC
    bulk_inorganicN[i,]<-bulk$inorganicN
    bulk_CO2[i,]<-bulk$CO2
    bulk_livingMicrobeN[i,]<-bulk$livingMicrobeN
    
  }
  
  ##This is a list of the information we want to get from this function
  Exported<-list(shared_inorganicN,litterbag,bulk,litterbag_uFastC,litterbag_uSlowC,litterbag_uNecroC,litterbag_uFastN,litterbag_uSlowN,litterbag_uNecroN,litterbag_livingMicrobeC,litterbag_inorganicN,litterbag_CO2, litterbag_livingMicrobeN,bulk_uFastC,bulk_uSlowC,bulk_uNecroC,bulk_livingMicrobeC,litterbag_pFastC,litterbag_pFastN,litterbag_pNecroC,litterbag_pNecroN,litterbag_pSlowC,litterbag_pSlowN)
  
  ##These are the names for each item in the list to be called later
  names(Exported)<-c("N inorg","litterbag","bulk","litterbag_uFastC","litterbag_uSlowC","litterbag_uNecroC","litterbag_uFastN","litterbag_uSlowN","litterbag_uNecroN","litterbag_livingMicrobeC","litterbag_inorganicN","litterbag_CO2", "litterbag_livingMicrobeN","bulk_uFastC","bulk_uSlowC","bulk_uNecroC","bulk_livingMicrobeC","litterbag_pFastC","litterbag_pFastN","litterbag_pNecroC","litterbag_pNecroN","litterbag_pSlowC","litterbag_pSlowN")
  
  return(Exported)
}



  
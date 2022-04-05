
#r Remove all functions, clear memory
rm(list=ls(all=TRUE)) 

#Load Packages
library(readr)

###Do you want to use new or old microbial params? Change to true if new params, false if OG params
NewParams=TRUE

##Load CORPSE functions
source("CORPSE_functions_variable.R")

##INPUT DATA
    #This reads in input variables for the function.
    CORPSE_full_spinup_bulk<-("Inputs/CORPSE_bulk.csv")
    litter_added<-("Inputs/CORPSE_litter.csv")
    soil_T<-paste("Inputs/CORPSE_soilT.csv",sep="")
    soil_VWC<-paste("Inputs/CORPSE_soilVWC.csv",sep="")
    
    ##Length of model run
    nyears<-1
    
    ##Defines how many litter types are being decomposed
    nspp<-4
    
    ##Running the model
    MainFunction<- CORPSE_MAIN(nyears,nspp,litter_added,CORPSE_full_spinup_bulk,soil_T,soil_VWC)


    

# Analysis and fast results ---------------------------------------------------

    #litter C remaining
    litter_C_remaining<-MainFunction$litterbag_uFastC+MainFunction$litterbag_uSlowC+MainFunction$litterbag_uNecroC+MainFunction$litterbag_livingMicrobeC
    
    #litter C protected
    litter_C_protected<-MainFunction$litterbag_pFastC+MainFunction$litterbag_pNecroC+MainFunction$litterbag_pSlowC
    
    #litter C unprotected
    litter_C_unprotected<-MainFunction$litterbag_uFastC+MainFunction$litterbag_uNecroC+MainFunction$litterbag_uSlowC
    
    #litter C respired
    litter_C_respired<-MainFunction$litterbag_CO2
    
   bulk_C_remaining<-MainFunction$bulk_uFastC+MainFunction$bulk_uSlowC+MainFunction$bulk_uNecroC+MainFunction$bulk_livingMicrobeC

##Results for termination of experiment at 12 weeks (day 84)
##div. by added 250mg of litter so output is % of added litter 
   print("respired")
   print(litter_C_respired[84,1:4]/.25)
   
   print("protected")
   print(litter_C_protected[84,1:4]/.25)

   print("unprotected")
   print(litter_C_unprotected[84,1:4]/.25)
   
  
   
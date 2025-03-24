################################
### Spatiotemporal analysis ####
################################

## load packages ---------------------------------------------------------------
library(INLA)
library(ggplot2)
library(cowplot)
library(ISOweek)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(forcats)


## read in data ----------------------------------------------------------------
# this dataset contains the estimated Rt for Malaysia at the state level between 2009-2019
# the variables matched to the data include school closure and climate variables over different lags
# please note that the cumulative incidence variable 'cum_inc' in the paper is not included in this data

data <- readRDS("data/Rt_variables.RDS") 


## load shapefile --------------------------------------------------------------
# source GADM: https://gadm.org/download_country.html

malaysia.1 <- readRDS("data/gadm36_MYS_1_sp.rds") 


## read in functions -----------------------------------------------------------

source("functions.R")


# run univariable analysis -----------------------------------------------------

all_variables <- colnames(data)[c(5,13:132)]
univar_results_table <- run_univar_analysis(data, variables = all_variables)

best_univariable_result<- best_univar(univar_results_table)
best_univariable_result

# plot univariable results -----------------------------------------------------

#choice must be all or only_significant
plot_univariable <- plot_univariable_beta(univar_results_table, choice = "all") 
plot_univariable


# run multivariable analysis ---------------------------------------------------
 
  # remove variables collinear with best variable from univariable analysis
  to_remove <- collinarity_check(data, all_variables, variable_chosen = best_univariable_result[[1]])
  step_variables <- all_variables[!all_variables %in% to_remove]
  to_remove <- step_variables[grep(substring(best_univariable_result[[1]], 1,7),step_variables)]
  step_variables <- step_variables[!step_variables %in% to_remove]
  

  # run step 1 of forward selection process
  step1 <- run_multivar_analysis(data, variables_to_try = step_variables, 
                               variables_chosen_already = best_univariable_result[[1]], stage = 1)
  
  # check waic decrease threshold reached
  (round(as.numeric(best_univariable_result[[2]])) - round(as.numeric(step1[[2]])))>4
  
  # remove variables collinear with best variable from step 1
  to_remove <- collinarity_check(data, all_variables, variable_chosen = step1[[1]])
  step_variables <- step_variables[!step_variables %in% to_remove]
  to_remove <- step_variables[grep(substring(step1[[1]], 1,7),step_variables)]
  step_variables <- step_variables[!step_variables %in% to_remove]
  
  
  # run step 2 of forward selection process
  step2 <- run_multivar_analysis(data, variables_to_try = step_variables, 
                                 variables_chosen_already = c(best_univariable_result[[1]], step1[[1]]), 
                                 stage = 2)
  
  # check waic decrease threshold reached
  (round(as.numeric(step1[[2]])) - round(as.numeric(step2[[2]])))>4
  
  # remove variables collinear with best variable from step 2
  to_remove <- collinarity_check(data, all_variables, variable_chosen = step2[[1]])
  step_variables <- step_variables[!step_variables %in% to_remove]
  to_remove <- step_variables[grep(substring(step2[[1]], 1,7),step_variables)]
  step_variables <- step_variables[!step_variables %in% to_remove]
  
  # run step 3 of forward selection process
  step3 <- run_multivar_analysis(data, variables_to_try = step_variables, 
                                 variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], step2[[1]]), 
                                 stage = 3)
  
  # check waic decrease threshold reached
  (round(as.numeric(step2[[2]])) - round(as.numeric(step3[[2]])))>4
  
  # remove variables collinear with best variable from step 3
  to_remove <- collinarity_check(data, all_variables, variable_chosen = step3[[1]])
  step_variables <- step_variables[!step_variables %in% to_remove]
  to_remove <- step_variables[grep(substring(step3[[1]], 1,7),step_variables)]
  step_variables <- step_variables[!step_variables %in% to_remove]
  
  # run step 4 of forward selection process
  step4 <- run_multivar_analysis(data, variables_to_try = step_variables, 
                                 variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], step2[[1]], step3[[1]]), 
                                 stage = 4)
  
  # check waic decrease threshold reached
  (round(as.numeric(step3[[2]])) - round(as.numeric(step4[[2]])))>4
  
  # remove variables collinear with best variable from step 4
  to_remove <- collinarity_check(data, all_variables, variable_chosen = step4[[1]])
  step_variables <- step_variables[!step_variables %in% to_remove]
  to_remove <- step_variables[grep(substring(step4[[1]], 1,7),step_variables)]
  step_variables <- step_variables[!step_variables %in% to_remove]
  
  # run step 5 of forward selection process
  step5 <- run_multivar_analysis(data, variables_to_try = step_variables, 
                                 variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], step2[[1]], step3[[1]], step4[[1]]), 
                                 stage = 5)
  
  # check waic decrease threshold reached
  (round(as.numeric(step4[[2]])) - round(as.numeric(step5[[2]])))>4
  
  # remove variables collinear with best variable from step 5
  to_remove <- collinarity_check(data, all_variables, variable_chosen = step5[[1]])
  step_variables <- step_variables[!step_variables %in% to_remove]
  to_remove <- step_variables[grep(substring(step5[[1]], 1,7),step_variables)]
  step_variables <- step_variables[!step_variables %in% to_remove]
  step_variables

  
# run final multivariable model and sample -------------------------------------
    model<- run_final_model(data) # this function has the variables pre coded
    waic <- model$waic$waic
    beta <- model$summary.fixed
    random <- model$summary.random
    
    
    samples <- sampling(model, data, nsamp=1000)
        
      
# plot final multivariable model results ---------------------------------------
plot_final_RE(random, choice="State_map", shapefile=malaysia.1) # choice must be State_line or State_map or Year
plot_final_RE(random, choice="State_line", shapefile=malaysia.1)
plot_final_RE(random, choice="Year", shapefile=malaysia.1)

plot_final_fit(samples, choice="Peninsular") # choice must be Peninsular or East
plot_final_fit(samples, choice="East")

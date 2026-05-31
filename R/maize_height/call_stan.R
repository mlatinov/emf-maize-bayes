#### Functions to Call the Stan models ####

#### M1 Plant Height ~ Treatment Model M1 ####
m1_stan_model <- function(model_data_maize_heigh, prior_only = 0){

  # Get the M1 Model 
  m1_model <- cmdstanr::cmdstan_model(stan_file = "Stan/maize_height/M1_maize_height.stan")

  # Sample 
  m1_sample <- m1_model$sample(
    data = list(
      N            = nrow(model_data_maize_heigh),
      treatment    = model_data_maize_heigh$treatment,
      plant_height = model_data_maize_heigh$plant_height,
      prior_only   = prior_only
    ),
    output_dir    = "stan_results/",
    iter_sampling = 1000,
    chains        = 4,
    seed          = 42    
  ) 
  return(m1_sample)
}

#### M2 Plant Height ~  Treatment × Day ####
m2_stan_model <- function(model_data_maize_heigh, prior_only = 0){

  # Get the model 
  m2_model <- cmdstanr::cmdstan_model(stan_file = "Stan/maize_height/M2_maize_height.stan")

  # Sample 
  m2_sample <- m2_model$sample(
    data = list(
      N = nrow(model_data_maize_heigh),
      treatment = model_data_maize_heigh$treatment,
      day       = model_data_maize_heigh$day,
      plant_height = model_data_maize_heigh$plant_height,
      prior_only   = prior_only
    ),
    output_dir    = "stan_results/",
    iter_sampling = 1000,
    chains        = 4,
    seed          = 42 
  )
  return(m2_sample)
}

#### M3 Plant Height ~ Pot x Day ####
m3_stan_model <- function(model_data_maize_heigh, prior_only = 0){

  # Get the model 
  m3_model <- cmdstanr::cmdstan_model(stan_file = "Stan/maize_height/M3_maize_height.stan")

  # Sample 
  m3_sample <- m3_model$sample(
    data = list(
      N   = nrow(model_data_maize_heigh),
      plant_height = model_data_maize_heigh$plant_height,
      day          = model_data_maize_heigh$day,
      pot_id       = model_data_maize_heigh$pot_id,
      prior_only = prior_only  
    ),
    output_dir    = "stan_results/",
    iter_sampling = 1000,
    chains        = 4,
    seed          = 42   
  )
  return(m3_sample)
}
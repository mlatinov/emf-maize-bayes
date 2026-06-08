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
      day       = model_data_maize_heigh$day_idx,
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
      day          = model_data_maize_heigh$day_idx,
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

#### M4 Gompertz per treatment ####
m4_stan_model <- function(model_data_maize_heigh, prior_only = 0, grid_points = 51){
  
  # Get the model 
  m4_model <- cmdstanr::cmdstan_model(stan_file = "Stan/maize_height/M4_maize_height.stan")

  # Sample 
  m4_sample <- m4_model$sample(
    data = list(
      N         = nrow(model_data_maize_heigh),
      treatment = model_data_maize_heigh$treatment,
      day_idx   = model_data_maize_heigh$day_idx,
      t         = model_data_maize_heigh$day,
      plant_height = model_data_maize_heigh$plant_height,
      G            = grid_points,
      prior_only   = prior_only 
    ),
    output_dir    = "stan_results/", 
    iter_sampling = 2000,
    chains        = 4,
    seed          = 42 
  )
  return(m4_sample)
}

#### Hierarchical cell-means, pots pooled within treatment ####
m5_stan_model <- function(
  model_data_maize_heigh,
  prior_only    = 0,
  iter_sampling = 4000,
  adapt_delta   = 0.95,
  chains        = 4
){
  # Get the model 
  model_5 <- cmdstanr::cmdstan_model(stan_file = "Stan/maize_height/M5_maize_height.stan")

  # Sample 
  m5_sample <- model_5$sample(
    data = list(
      N  = nrow(model_data_maize_heigh),
      plant_height  = model_data_maize_heigh$plant_height,
      treatment     = model_data_maize_heigh$treatment,
      day_idx       = model_data_maize_heigh$day_idx,
      pot_idx       = model_data_maize_heigh$pot_id,
      pot_treatment = c(rep(1, 6), rep(2, 4), rep(3, 6)),
      prior_only    = prior_only   
    ),
    output_dir    = "stan_results/",  
    iter_sampling = iter_sampling,
    adapt_delta   = adapt_delta,
    chains        = chains,
    seed = 42
  )
  return(m5_sample)
}

#### Hierarchical Gompertz, per-pot curves pooled within treatment ####
m6_stan_model <- function(
  model_data_maize_heigh, 
  prior_only    = 0,
  iter_sampling = 4000,
  chains        = 4,
  adapt_delta   = 0.95,
  grid_points   = 51
){

  # Get the model 
  m6_model <- cmdstanr::cmdstan_model(stan_file = "Stan/maize_height/M6_maize_height.stan")

  # Sample 
  m6_sample <- m6_model$sample(
    data = list(
      N  = nrow(model_data_maize_heigh),
      plant_height  = model_data_maize_heigh$plant_height,
      day_idx       = model_data_maize_heigh$day_idx,
      t             = model_data_maize_heigh$day,
      pot_idx       = model_data_maize_heigh$pot_id,
      pot_treatment = c(rep(1, 6), rep(2, 4), rep(3, 6)),
      prior_only    = prior_only,
      G             = grid_points
    ),
    output_dir    = "stan_results/",
    iter_sampling = iter_sampling,
    chains        = chains,
    adapt_delta   = adapt_delta,
    seed = 42
  )
  return(m6_sample)
}

#### Function to Construct effect size curves 
get_effect_size_grid <- function(var, g){
  var_log_range <- range(log(var))
  var_grid      <- seq(
    (var_log_range[1] - mean(log(var))) / sd(log(var)),
    (var_log_range[2] - mean(log(var))) / sd(log(var)),
  length.out = g
  )
  return(var_grid)
}

#### Function to Call Stan model M1 Maize BB ####
bb_model <- function(
  data_bb, 
  prior_only   = 0, 
  stan_file    = "Stan/maize_biochem_biomass/M1_maize_bb.stan", 
  iter         = 2000, 
  adapt_delta  = 0.95,
  grid_points  = 50
){

  ## Prepare the data for Stan 
  data             <- data_bb$data
  pot_treatment_id <- data_bb$additional$pot_treatment_id

  ## Construct effect size curves
  sod_z_grid  <- get_effect_size_grid(data_bb$data$sod_dw,   g = grid_points)
  cat_z_grid  <- get_effect_size_grid(data_bb$data$cat_dw,   g = grid_points)
  h2o2_z_grid <- get_effect_size_grid(data_bb$data$h2_o2_dw, g = grid_points)
  
  ## Combine everything in Stan data list to pass to cmdstanr
  stan_data <- list(
    N = nrow(data),
    G = grid_points,
    sod_z_grid  = sod_z_grid,
    h202_z_grid = cat_z_grid,
    cat_z_grid  = h2o2_z_grid,
    prior_only = prior_only,
    pot_treatment_id = pot_treatment_id,
    pot_id           = data$pot_id,
    treatment_id = data$treatment_id,
    day_id       = data$day_id,
    cell_id      = data$cell_id,
    sod_dw     = data$sod_dw,
    cat_dw     = data$cat_dw,
    sugar_dw   = data$reduced_sugar_dw,
    trolox_dw  = data$trolox_umol_dw,
    h2o2_dw    = data$h2_o2_dw,
    mda_dw     = data$mda_dw, 
    water_content_dw = data$water_content_dw
  )

  ## Compile the stan model 
  stan_model <- cmdstanr::cmdstan_model(stan_file = stan_file)

  ## Sample from the model 
  sample <- stan_model$sample(
    data = stan_data,
    iter_sampling = iter,
    chains        = 4,
    output_dir    = "stan_results/",
    seed          = 42,
    adapt_delta   = adapt_delta
  )
  return(list(sample = sample, stan_data = stan_data))
}
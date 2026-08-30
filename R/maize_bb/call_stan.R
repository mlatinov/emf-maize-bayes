#### Function to Call Stan model M1 Maize BB ####
bb_model <- function(
  data_bb, 
  prior_only   = 0, 
  stan_file    = "Stan/maize_biochem_biomass/M1_maize_bb.stan", 
  iter         = 2000, 
  adapt_delta  = 0.95
){

  ## Prepare the data for Stan 
  data             <- data_bb$data
  pot_treatment_id <- data_bb$additional$pot_treatment_id

  stan_data <- list(
    N = nrow(data),
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
  return(sample)
}
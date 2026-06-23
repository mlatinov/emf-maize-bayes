#### Function to Call Stan Model for JIP analysis ####

#### GammaRC Model 1 Cell-Means Model Without Heirarcle or Pooling 
m1_gamma_rc_model <- function(model_data_jip, prior_only = 0){

  # Compile Stan model 
  m1_gamma_rc <- cmdstanr::cmdstan_model(stan_file = "Stan/maize_jip/M1_maize_jip_gammarc.stan")

  # Sample from the model 
  m1_gamma_rc_sample <- m1_gamma_rc$sample(
    data = list(
      N = nrow(model_data_jip),
      prior_only = prior_only,
      gammaRC    = model_data_jip$gammaRC,
      day_idx    = model_data_jip$day_idx,
      treatment  = model_data_jip$treatment
    ),
    output_dir    = "stan_results/",
    chains        = 4,
    iter_sampling = 2000,
    seed          = 42 
  )
  # Return the mdodel Output 
  return(m1_gamma_rc_sample)
}

#### GammaRC2 Model Hierarcle Model with Parial Pooling ####
m2_gamma_rc_model <- function(model_data_jip, prior_only = 0) {
  
  ## Compile the model 
  m2_gamma_rc <- cmdstanr::cmdstan_model(stan_file = "Stan/maize_jip/M2_maize_jip_gammarc.stan")

  ## Initial Paramter Function to Start 
  init_fn <- function() list(
    eta     = matrix(-1, 3, 3), # logit of ~0.27
    kappa_d = rep(50, 3),       # concentration
    tau_t   = rep(0.1, 3),      # small pot SD
    zp      = rep(0, 16)        # pots start at cell mean
  )
  
  # Sample from the model 
  sample_gamma_rc <- m2_gamma_rc$sample(
    data = list(
      N = nrow(model_data_jip),
      prior_only = prior_only,
      gammaRC    = model_data_jip$gammaRC,
      day_idx    = model_data_jip$day_idx,
      treatment  = model_data_jip$treatment,
      pot_idx    = model_data_jip$pot_idx,
      pot_treatment = c(rep(1, 6), rep(2, 4), rep(3, 6))
    ),
    output_dir    = "stan_results/",
    iter_sampling = 3000,
    chains        = 4,
    adapt_delta   = 0.95,
    init = init_fn,
    seed = 123
  )
  # Return the model object 
  return(sample_gamma_rc)
}
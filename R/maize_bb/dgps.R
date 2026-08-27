#### Function to simulate primaries LogNormal ####
## This function uses the same mechanism with to simulate sod_dw, cat_dw, reduced_sugar_dw and trolox_umol_dw
## This can be achieved by changing the cell mean ands sd as well as zeta rate 
simulate_primaries <- function(
  n_treatment     = 3,
  n_day           = 3,
  sdlog_obs       = 0.08,     # residual SD on the LOG scale
  cell_mean       = 220,      # mean SOD activity
  cell_sd_log     = 0.20,     # spread of cell means
  zeta_rate       = 15,
  primary_name    = "sod_dw"       
){

  ## Simulate the experimental design settings ##
  id <- rsims::make_design(
    nested  = list(treatment = c(Control = 6, Sham = 4, EMF = 6)),
    crossed = list(day = c(14, 21, 28))
  )
  n_pot <- length(unique(id$unit_id))

  ## Model parameters, all on the LOG scale ##
  cell_matrix <- matrix(
    data = rnorm(n_treatment * n_day, mean = log(cell_mean), sd = cell_sd_log),
    nrow = n_treatment,
    ncol = n_day
  )

  zeta_t <- rexp(n_treatment, rate = zeta_rate)   # treatment-specific pot-scatter SD, LOG scale
  z_pot  <- rnorm(n_pot, mean = 0, sd = 1)        # raw unit-normal deviate

  ## Lookup:
  pot_treatment_id <- id$treatment_idx[!duplicated(id$unit_id)]

  ## Non-centered pot offset: SD * raw deviate, on the log scale
  u_pot <- zeta_t[pot_treatment_id] * z_pot

  ## Linear predictor 
  meanlog <- cell_matrix[id$cell_idx] + u_pot[id$unit_id]

  ## Sample from LogNormal
  primary <- rlnorm(nrow(id), meanlog = meanlog, sdlog = sdlog_obs) 

  ## Return both the simulated data AND the ground truth
  list(
    data = data.frame(id, primary) %>% rename(!!primary_name := primary),
    truth = list(
      cell_matrix   = cell_matrix,   # log-scale cell means
      zeta_t        = zeta_t,
      z_pot         = z_pot,
      u_pot         = u_pot,
      meanlog       = meanlog        
    )
  )
}

#### Function to simulate secondaries LogNormal ####
## This Simulates h2_o2 it uses the simulate_primaries for sod_dw, cat_dw, trolox_umol_dw 
simulate_h2_02 <- function(
  beta_sod  = 0.5,
  beta_cad  = 0.05,
  beta_trolox = 0.08,
  mean_cell    = 110,
  mean_cell_sd = 0.01,
  zeta_rate    = 50,
  sdlog_obs    = 0.02,
  seed         = 42
){

  ## Simulate the experimental design settings ##
  id <- rsims::make_design(
    nested  = list(treatment = c(Control = 6, Sham = 4, EMF = 6)),
    crossed = list(day = c(14, 21, 28))
  )
  n_pot            <- length(unique(id$unit_id))
  pot_treatment_id <- id$treatment_idx[!duplicated(id$unit_id)]

  ## Simulate the necessary primaries ## 
  set.seed(seed)
  # SOD 
  sod_dw <- simulate_primaries(  
    sdlog_obs       = 0.08,     
    cell_mean       = 180,      
    cell_sd_log     = 0.25,     
    zeta_rate       = 20,
    primary_name    = "sod_dw" 
  )

  # CAT 
  cat_dw <- simulate_primaries(  
    sdlog_obs       = 0.08,     
    cell_mean       = 2,      
    cell_sd_log     = 0.20,     
    zeta_rate       = 10,
    primary_name    = "cat_dw" 
  )

  # TEAC 
  trolox_umol_dw <- simulate_primaries(
    sdlog_obs       = 0.08,     
    cell_mean       = 130,      
    cell_sd_log     = 0.25,     
    zeta_rate       = 15,
    primary_name    = "trolox_umol_dw" 
  )

  ## Model Parameters 
  cell_matrix <- matrix(data = rnorm(9, mean = log(mean_cell), sd = mean_cell_sd), 3, 3)
  zeta_t      <- rexp(3, rate = zeta_rate)   # treatment-specific pot-scatter SD, LOG scale
  z_pot       <- rnorm(16, mean = 0, sd = 1) # raw unit-normal deviate

  ## Non-centered pot offset: SD * raw deviate, on the log scale
  u_pot <- zeta_t[pot_treatment_id] * z_pot

  ## Standardize the simulated primaries 
  mu_sod_z <- scale(sod_dw$truth$meanlog)
  mu_cat_z <- scale(cat_dw$truth$meanlog)
  mu_trolox_z <- scale(trolox_umol_dw$truth$meanlog)

  ## Linear predictor 
  meanlog <- (
    cell_matrix[id$cell_idx] + u_pot[id$unit_id]
    + beta_sod * mu_sod_z
    + beta_cad * mu_cat_z
    + beta_trolox * mu_trolox_z
  )

  ## Sample from LogNormal
  h2_o2 <- rlnorm(nrow(id), meanlog = meanlog, sdlog = sdlog_obs) 

  ## Return the H2o2 the design and the primaries + truth
  sim_data <- data.frame(
    id,
    cat_dw = cat_dw$data$cat_dw,
    trolox_umol_dw = trolox_umol_dw$data$trolox_umol_dw,
    sod_dw         = sod_dw$data$sod_dw,
    h2_o2_dw       = h2_o2
  )
  truth <- list(
    meanlog     = meanlog,
    cell_matrix = cell_matrix,
    zeta_t      = zeta_t,
    z_pot       = z_pot
  )
  primaries <- list(sod_dw, cat_dw, trolox_umol_dw)

  return(list(
    sim = sim_data, 
    truth = truth, 
    primaries = primaries)
  )
}

#### Function to simulate MDA ###
simulate_mda <- function(
  cell_sd   = 0.02,
  cell_mean = 300,
  sdlog_obs = 0.01,
  zeta_rate = 30,
  beta_h2o2 = 0.5,
  partial_beta_h2o2 = 0.5,
  sd_partial_beta_h2o2 = 0.2,
  partial_pooling_h2o2 = FALSE
){
  
  ## Simulate the experimental design settings ##
  id <- rsims::make_design(
    nested  = list(treatment = c(Control = 6, Sham = 4, EMF = 6)),
    crossed = list(day = c(14, 21, 28))
  )
  n_pot            <- length(unique(id$unit_id))
  pot_treatment_id <- id$treatment_idx[!duplicated(id$unit_id)]

  ## Simulate the h2o2 ## 
  h2_o2_data <- simulate_h2_02() 
  mu_h2o2_z  <- scale(h2_o2_data$truth$meanlog)

  ## Model parameters
  cell_matrix <- matrix(data = rnorm(9, mean = log(cell_mean), sd = cell_sd))
  zeta_t      <- rexp(3, rate = zeta_rate)   # treatment-specific pot-scatter SD, LOG scale
  z_pot       <- rnorm(16, mean = 0, sd = 1) # raw unit-normal deviate
  u_pot       <- zeta_t[pot_treatment_id] * z_pot
  
  # Allow for Partial Pooling of the beta_h2o2 #
  if(partial_pooling_h2o2){
    # Pool beta h202
    random_beta_h2o2 <- partial_beta_h2o2 + sd_partial_beta_h2o2 * rnorm(3, mean = 0, sd = 1)

    # Linear Predictor
    mean_log <- (
      cell_matrix[id$cell_idx] + u_pot[id$unit_id] 
      + mu_h2o2_z * random_beta_h2o2[id$treatment_idx] 
    )
  }else{

    # Linear Predictor
    mean_log <- (
      cell_matrix[id$cell_idx] + u_pot[id$unit_id] 
      + mu_h2o2_z * beta_h2o2
    )
  }

  # Sample from Lognormal Distribution 
  mda_dw <- rlnorm(nrow(id), meanlog = mean_log, sdlog = sdlog_obs)

  # Combine and return list of the data and truth
  mda_data <- data.frame(id, mda_dw)
  
  sim_data <- data.frame(
    id,
    mda_dw = mda_dw,
    sod_dw = h2_o2_data$sim$sod_dw,
    cat_dw = h2_o2_data$sim$cat_dw,
    trolox_umol_dw = h2_o2_data$sim$trolox_umol_dw,
    h2_o2_dw       = h2_o2_data$sim$h2_o2_dw  
  )

  truth <- list(
    meanlog     = meanlog,
    cell_matrix = cell_matrix,
    zeta_t      = zeta_t,
    z_pot       = z_pot
  )

  return(list(
    sim_data   = sim_data,
    h2_o2_data = h2_o2_data,
    truth      = truth
  ))
}
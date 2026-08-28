#### Function to simulate primaries LogNormal ####
## This function uses the same mechanism with to simulate sod_dw, cat_dw, reduced_sugar_dw and trolox_umol_dw
## This can be achieved by changing the cell mean ands sd as well as zeta rate 
## Note : Currently All the functions in this file simulates DW covariates which is the primary objective ##
simulate_primaries <- function(
  n_treatment     = 3,
  n_day           = 3,
  lambda_sdlog_obs = 10.5,    # residual SD on the LOG scale
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
  id$pot_id <- id$unit_id
  n_pot <- length(unique(id$unit_id))

  ## Model parameters, all on the LOG scale ##
  cell_matrix <- matrix(
    data = rnorm(n_treatment * n_day, mean = log(cell_mean), sd = cell_sd_log),
    nrow = n_treatment,
    ncol = n_day
  )

  zeta_t <- rexp(n_treatment, rate = zeta_rate)   # treatment-specific pot-scatter SD, LOG scale
  z_pot  <- rnorm(n_pot, mean = 0, sd = 1)        # raw unit-normal deviate
  sdlog_obs <- rexp(3, rate = lambda_sdlog_obs)

  ## Lookup:
  pot_treatment_id <- id$treatment_idx[!duplicated(id$unit_id)]

  ## Non-centered pot offset: SD * raw deviate, on the log scale
  u_pot <- zeta_t[pot_treatment_id] * z_pot

  ## Linear predictor 
  meanlog <- cell_matrix[id$cell_idx] + u_pot[id$unit_id]

  ## Sample from LogNormal
  primary <- rlnorm(nrow(id), meanlog = meanlog, sdlog = sdlog_obs[id$day_idx]) 

  ## Return both the simulated data AND the ground truth
  list(
    data = data.frame(id, primary) %>% rename(!!primary_name := primary),
    truth = list(
      cell_matrix   = cell_matrix,   # log-scale cell means
      zeta_t        = zeta_t,
      z_pot         = z_pot,
      u_pot         = u_pot,
      meanlog       = meanlog,
      sdlog_obs     = sdlog_obs
    ),
    additional = list(pot_treatment_id = pot_treatment_id)    
  )
}

#### Function to simulate H202 DW uses the simulate_primaries for sod_dw, cat_dw, trolox_umol_dw #### 
simulate_h2_02 <- function(
  beta_sod    = 0.21,
  beta_cad    = 0.05,
  beta_trolox = 0.05,
  lambda_sdlog_obs = 5.7,    
  mean_cell   = 110,       
  mean_cell_sd = 0.01,
  zeta_rate    = 50,
  seed         = 42
){
  ## Simulate the experimental design settings ##
  id <- rsims::make_design(
    nested  = list(treatment = c(Control = 6, Sham = 4, EMF = 6)),
    crossed = list(day = c(14, 21, 28))
  )
  id$pot_id <- id$unit_id
  n_pot            <- length(unique(id$unit_id))
  pot_treatment_id <- id$treatment_idx[!duplicated(id$unit_id)]

  ## Simulate the necessary primaries ## 
  set.seed(seed)

  # SOD 
  sod_dw <- simulate_primaries(    
    cell_mean       = 180,      
    cell_sd_log     = 0.25,     
    zeta_rate       = 20,
    lambda_sdlog_obs = 10.5,
    primary_name    = "sod_dw" 
  )

  # CAT 
  cat_dw <- simulate_primaries(  
    cell_mean       = 2,      
    cell_sd_log     = 0.20,     
    zeta_rate       = 10,
    lambda_sdlog_obs = 3.4,
    primary_name    = "cat_dw" 
  )

  # TEAC 
  trolox_umol_dw <- simulate_primaries(
    cell_mean       = 130,      
    cell_sd_log     = 0.25,     
    zeta_rate       = 15,
    lambda_sdlog_obs = 8.3,
    primary_name     = "trolox_umol_dw" 
  )

  ## Model Parameters 
  cell_matrix <- matrix(data = rnorm(9, mean = log(mean_cell), sd = mean_cell_sd), 3, 3)
  zeta_t      <- rexp(3, rate = zeta_rate)   # treatment-specific pot-scatter SD, LOG scale
  z_pot       <- rnorm(16, mean = 0, sd = 1) # raw unit-normal deviate
  sdlog_obs   <- rexp(3, rate = lambda_sdlog_obs) 
  
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
  h2_o2 <- rlnorm(nrow(id), meanlog = meanlog, sdlog = sdlog_obs[id$day_idx]) 

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
    z_pot       = z_pot,
    sdlog_obs   = sdlog_obs
  )
  primaries <- list(sod_dw, cat_dw, trolox_umol_dw)

  return(list(
    data  = sim_data, 
    truth = truth, 
    primaries = primaries,
    additional = list(pot_treatment_id = pot_treatment_id)
  ))
}

#### Function to simulate MDA DW ###
simulate_mda <- function(
  beta_h2o2  = 0.18,
  lambda_sdlog_obs = 5.3,     
  cell_sd    = 0.02,
  cell_mean  = 300,
  zeta_rate  = 30,
  partial_beta_h2o2 = 0.5,
  sd_partial_beta_h2o2 = 0.2,
  partial_pooling_h2o2 = FALSE
){
  
  ## Simulate the experimental design settings ##
  id <- rsims::make_design(
    nested  = list(treatment = c(Control = 6, Sham = 4, EMF = 6)),
    crossed = list(day = c(14, 21, 28))
  )
  id$pot_id        <- id$unit_id
  n_pot            <- length(unique(id$unit_id))
  pot_treatment_id <- id$treatment_idx[!duplicated(id$unit_id)]

  ## Simulate the h2o2 ## 
  h2_o2_data <- simulate_h2_02() 
  mu_h2o2_z  <- scale(h2_o2_data$truth$meanlog)

  ## Model parameters
  cell_matrix <- matrix(data = rnorm(9, mean = log(cell_mean), sd = cell_sd),3, 3)
  zeta_t      <- rexp(3, rate = zeta_rate)   # treatment-specific pot-scatter SD, LOG scale
  z_pot       <- rnorm(16, mean = 0, sd = 1) # raw unit-normal deviate
  u_pot       <- zeta_t[pot_treatment_id] * z_pot
  sdlog_obs   <- rexp(3, rate = lambda_sdlog_obs)

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
  mda_dw <- rlnorm(nrow(id), meanlog = mean_log, sdlog = sdlog_obs[id$day_idx])

  # Combine and return list of the data and truth
  sim_data <- data.frame(
    id,
    mda_dw = mda_dw,
    sod_dw = h2_o2_data$data$sod_dw,
    cat_dw = h2_o2_data$data$cat_dw,
    trolox_umol_dw = h2_o2_data$data$trolox_umol_dw,
    h2_o2_dw       = h2_o2_data$data$h2_o2_dw  
  )
  truth <- list(
    meanlog     = mean_log,
    cell_matrix = cell_matrix,
    zeta_t      = zeta_t,
    z_pot       = z_pot,
    sdlog_obs   = sdlog_obs
  )
  return(list(
    data       = sim_data,
    h2_o2_data = h2_o2_data,
    truth      = truth,
    additional = list(pot_treatment_id = pot_treatment_id)    
  ))
}

#### Function to Simulate Biomass (Leaf Water Content Relative to DW) From MDA DW and Resid Sugars DW ####
simulate_biomass <- function(
  cell_mean  = 9.5,
  cell_sd    = 0.05,
  lambda_sdlog_obs  = 11,
  zeta_rate  = 15,
  beta_sugar = 0.05,
  beta_mda   = 0.01,
  seed       = 42
){
  ## Simulate the experimental design settings ##
  set.seed(seed)
  id <- rsims::make_design(
    nested  = list(treatment = c(Control = 6, Sham = 4, EMF = 6)),
    crossed = list(day = c(14, 21, 28))
  )
  id$pot_id        <- id$unit_id
  n_pot            <- length(unique(id$unit_id))
  pot_treatment_id <- id$treatment_idx[!duplicated(id$unit_id)]

  ## Simulate MDA 
  mda_data <- simulate_mda()

  ## Simulate Resid Sugars as a primary
  sugar_data <- simulate_primaries(
    lambda_sdlog_obs = 3.5,
    cell_mean   = 80,
    cell_sd_log = 0.02,
    zeta_rate   = 15,
    primary_name = "reduced_sugar_dw"
  ) 

  # Scale the meanlogs
  mu_mda_z   <- scale(mda_data$truth$meanlog)
  mu_sugar_z <- scale(sugar_data$truth$meanlog)

  ## Biomass Model parameters
  cell_matrix <- matrix(data = rnorm(9, mean = log(cell_mean), sd = cell_sd), 3, 3)
  zeta_t      <- rexp(3, rate = zeta_rate)   
  z_pot       <- rnorm(16, mean = 0, sd = 1) 
  u_pot       <- zeta_t[pot_treatment_id] * z_pot
  sdlog_obs   <- rexp(3, rate = lambda_sdlog_obs) 

  ## Linear Predictor
  meanlog <- (
    cell_matrix[id$cell_idx] + u_pot[id$unit_id] 
    + mu_mda_z   * beta_mda
    + mu_sugar_z * beta_sugar
  )

  ## Sample from Lognormal Distribution 
  water_content_dw <- rlnorm(nrow(id), meanlog = meanlog, sdlog = sdlog_obs[id$day_idx])

  ## Combine and Return the data and truth
  sim_data <- data.frame(
    id,
    water_content_dw = water_content_dw,
    mda_dw           = mda_data$data$mda_dw,
    sod_dw           = mda_data$data$sod_dw,
    cat_dw           = mda_data$data$cat_dw,
    trolox_umol_dw = mda_data$data$trolox_umol_dw,
    h2_o2_dw       = mda_data$data$h2_o2_dw,
    reduced_sugar_dw = sugar_data$data$reduced_sugar_dw
  )
  truth <- list(
    meanlog     = meanlog,
    cell_matrix = cell_matrix,
    zeta_t      = zeta_t,
    z_pot       = z_pot,
    sdlog_obs   = sdlog_obs
  )
  primary <- list(
      sugar_data = sugar_data,
      mda_data   = mda_data
  )
  
  return(list(
    data     = sim_data,
    truth    = truth,
    primary  = primary,
    additional = list(pot_treatment_id = pot_treatment_id)    
  ))
}
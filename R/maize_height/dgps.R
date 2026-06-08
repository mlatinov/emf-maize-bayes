
#### Functions for Data Generative Processes for Maize Plant Height ####

#### Simple Gaussian Model Generation Process ####
maize_h_simple_gaussian_dgp <- function(
  n        = 700,
  meanlog  = c(log(15), log(15), log(18)),  
  sdlog    = 0.2,                            
  prob     = c(6, 4, 6) / 16  
){

  # Simulate the treatment 
  treatment <- sample(1:3, size = n, replace = TRUE, prob = prob)

  # Generate outcome from LogNormal Distribution 
  height  <- rlnorm(n, meanlog = meanlog[treatment], sdlog = sdlog)

  # Return a tibble 
  return(tibble(treatment = treatment, height = height))
} 

#### M2 — Plant Height ~ Treatment × Day ####
maize_h_m2_dgp <- function(
  n = 700,
  meanlog_cell = matrix(
    # Day 7,       14,       20,     27
    c(log(10), log(20), log(30), log(40),   # Control
      log(8),  log(18), log(26), log(35),   # Sham 
      log(11), log(22), log(33), log(42)    # EMF 
    ),  
    nrow = 3, ncol = 4, byrow = TRUE
  ),
  sdlog          = 0.2,
  prob_treatment = c(6, 4, 6) / 16,
  prob_day       = c(255, 249, 169, 95) / 768
){
  # Simulate treatment and Day
  treatment <- sample(1:3, size = n, replace = TRUE, prob = prob_treatment)
  day       <- sample(1:4, size = n, replace = TRUE, prob = prob_day)

  # Generate Plant Heights 
  meanlog   <- meanlog_cell[cbind(treatment, day)]  
  height    <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)

  # Return a Tibble 
  return(tibble(treatment = treatment, day = day, height = height))
}

#### M3 Plant Height ~ Pot x Day ####
maize_h_m3_dgp <- function(
  n = 700,
  mu_meanlog  = 3,
  sd_mean_log = 0.3, 
  prob_day = c(255, 249, 169, 95) / 768,
  prob_pot = c(53,49,44,52,41,43,54,51,45,54,55,52,40,45,41,49 ) / 768,
  sdlog = 0.2
){
  # Simulate Day and Pot
  pot <- sample(1:16, size = n, replace = TRUE, prob = prob_pot)
  day <- sample(1:4, size = n, replace = TRUE, prob = prob_day)

  # Mean Log Matrix 
  meanlog_cell <- matrix(
    rnorm(16 * 4, mean = mu_meanlog, sd = sd_mean_log),
    nrow = 16,
    ncol  = 4
  )
  
  # Generate Plant Heights 
  meanlog <- meanlog_cell[cbind(pot, day)]
  height  <- rlnorm(n, meanlog = meanlog, sdlog = sdlog) 

  # Return a Tibble 
  return(tibble(pot_id = pot, day = day, height = height))
}

#### M4 Gompertz per treatment, no pots ####
maize_h_m4_dgp <- function(
  n = 700,
  prob_treatment = c(6, 4, 6) / 16,
  prob_day       = c(255, 249, 169, 95) / 768
){

  # Treatment + day index
  treatment <- sample(1:3, n, replace = TRUE, prob = prob_treatment)
  day       <- sample(1:4, n, replace = TRUE, prob = prob_day)

  # Real time
  time_values <- c(7, 14, 20, 27)
  t <- time_values[day]

  # Parameters
  A_t   <- rlnorm(n = 3, meanlog =  log(40),  sdlog = 0.3)
  k_t   <- rlnorm(n = 3, meanlog =  log(0.15),sdlog =  0.5)
  tao_t <- rnorm(n = 3, mean = 9, sd = 3)
  sigma_day <- rexp(n = 4, rate = 2)

  # Gompertz Equation for the mean 
  mu_i <- A_t[treatment] * exp(-exp(-k_t[treatment] * (t - tao_t[treatment])))

  # Simulate Plant height 
  plant_height <- rlnorm(n, log(mu_i), sigma_day[day])

  return(
    data.frame(
      treatment = treatment,
      day_index = day,
      day = t,
      height = plant_height
    )
  )
}

#### Hierarchical cell-means, pots pooled within treatment ####
maize_h_m5_dgp <- function(
  n = 700,
  prob_treatment = c(6, 4, 6) / 16,
  prob_day       = c(255, 249, 169, 95) / 768,
  prob_pot       = c(53,49,44,52,41,43,54,51,45,54,55,52,40,45,41,49 ) / 768
){

  # Simulate Treatment Day Index and Pot Index 
  pot_treatment <- c(rep(1, 6),rep(2, 4),rep(3, 6))
  pot_idx   <- sample(1:16, size = n, replace = TRUE, prob = prob_pot)
  treatment <- pot_treatment[pot_idx] 
  day_idx   <- sample(1:4, size = n, replace = TRUE, prob = prob_day)

  # Residual Scale One per day 
  sigma_d <- rexp(n = 4, rate = 2)

  # Pot Scatter One per treatment 
  zeta_t  <- rexp(n = 3 ,rate =  10)

  # Population Level Means Matrix by treatment and day 
  alpha_dt <- matrix(
    rnorm(n = 12, mean = 3, sd = 1),
    nrow = 3,
    ncol = 4
  )

  # Unit-scale pot deviation
  zeta_p_bar <- rnorm(n = 16, mean = 0, sd = 1) 

  # Pot-level offsets
  zp <- zeta_t[pot_treatment] * zeta_p_bar

  # Observation-level likelihood
  mu     <- alpha_dt[cbind(treatment, day_idx)] + zp[pot_idx]
  height <- rlnorm(n = n, meanlog = mu, sdlog = sigma_d[day_idx])

  # Return the data and the paramters 
  gen_data <- tibble(
    treatment = treatment,
    day_idx   = day_idx,
    pot_id    = pot_idx,
    plant_height = height
  )

  # Paramtemers 
  params <- list(
    sigma_d  = sigma_d,
    zeta_t   = zeta_t,
    alpha_dt = alpha_dt,
    zp = zp,
    mu = mu
  )

  return(list(
    data = gen_data,
    parameters = params
  ))
}

#### Hierarchical Gompertz, per-pot curves pooled within treatment ####
maize_h_m6_dgp <- function(
  n              = 768,
  pot_treatment  = c(rep(1, 6), rep(2, 4), rep(3, 6)),
  prob_pot       = c(53, 49, 44, 52, 41, 43, 54, 51, 45, 54, 55, 52, 40, 45, 41, 49) / 768,
  prob_day       = c(255, 249, 169, 95) / 768,
  
  # Population curve parameters per treatment 
  a_t       = c(log(39), log(34), log(37)),        # log-asymptote: Control, Sham, EMF
  b_t       = c(log(0.15), log(0.18), log(0.18)),  # log-rate
  c_t       = c(8.9, 8.9, 8.0),                    # inflection time
  
  # Pot-scatter SDs 
  xi_A_t    = c(0.05, 0.06, 0.05),
  xi_K_t    = c(0.10, 0.10, 0.10),
  xi_tau_t  = c(0.1, 0.1, 0.1),
  
  # Per-day residual SD 
  sigma_d   = c(0.47, 0.34, 0.32, 0.24)
){
  # Indices 
  pot_idx     <- sample(1:16, size = n, replace = TRUE, prob = prob_pot)
  day_idx     <- sample(1:4, size = n, replace = TRUE, prob = prob_day)
  treatment   <- pot_treatment[pot_idx]
  time_values <- c(7, 14, 20, 27)
  day         <- time_values[day_idx]
  
  # Unit-scale pot deviations 
  z_A   <- rnorm(n = 16, mean = 0, sd = 1)
  z_K   <- rnorm(n = 16, mean = 0, sd = 1)
  z_tau <- rnorm(n = 16, mean = 0, sd = 1)
  
  # Pot-level curve parameters 
  log_A_p <- a_t[pot_treatment] + xi_A_t[pot_treatment]   * z_A
  log_K_p <- b_t[pot_treatment] + xi_K_t[pot_treatment]   * z_K
  tau_p   <- c_t[pot_treatment] + xi_tau_t[pot_treatment] * z_tau
  
  # Convert to natural scale for the Gompertz curve
  A_p <- exp(log_A_p)
  K_p <- exp(log_K_p)
  
  # Gompertz curve
  mu_i <- A_p[pot_idx] * exp(-exp(-K_p[pot_idx] * (day - tau_p[pot_idx])))
  
  # LogNormal observation 
  height <- rlnorm(n = n, meanlog = log(mu_i), sdlog = sigma_d[day_idx])
  
  # Return data
  return(list(
    data = tibble::tibble(
      day_idx      = day_idx,
      day          = day,
      pot_id       = pot_idx,
      treatment    = treatment,
      plant_height = height
    ),
    params = list(
      a_t = a_t, b_t = b_t, c_t = c_t,
      xi_A_t = xi_A_t, xi_K_t = xi_K_t, xi_tau_t = xi_tau_t,
      sigma_d = sigma_d,
      log_A_p = log_A_p, log_K_p = log_K_p, tau_p = tau_p
    )
  ))
}



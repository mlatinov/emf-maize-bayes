
#### Function to for Simulating Data Generative Processes to test the model ####
set.seed(42)

#### Cell- Means model Without Pooling ####
dgp_gammaRC_jip_1 <- function(
  n = 188,
  prob_treatment = c(68, 72, 48) / 188,
  prob_day       = c(64, 64, 60) / 188
){
  # Sumulate Day and  Treatment
  day       <- sample(1:3, size = n, replace = TRUE, prob = prob_day)
  treatment <- sample(1:3, size = n, replace = TRUE, prob = prob_treatment)

  # Simulate Prior Matrix from Normal Distribution for each d and t 
  eta <- matrix(rnorm(n = 9, mean = -1, sd = 1),nrow = 3,ncol = 3,dimnames = list(treatment = 1:3, day = 1:3))
  mu_mat <- plogis(eta)
  
  # For every day D allow seperate kappa variation 
  kappa_d <- rexp(n = 3, rate = 0.01)

  # Per obervation pair
  mu_i    <- mu_mat[cbind(treatment, day)]
  kappa_i <- kappa_d[day]
  
  # Sample from Beta Distribution the Gamma RC 
  gammaRC <- rbeta(
    n = n, 
    shape1 = mu_i * kappa_i,
    shape2 = (1 - mu_i) * kappa_i
  )
  # Return the simulated Data and the Paramters 
  return(list(
    data = tibble(
      day_idx   = day,
      treatment = treatment,
      gammaRC   = gammaRC
    ),
    parameters = tibble(
      day_idx = 1:3,
      kappa_per_day = kappa_d,
      cell_means_logit = eta, 
      cell_means_prob  = mu_mat
    )
  ))
}

#### Cell-Means Model With Pot Pooling ####
dgp_gammaRC_jip_2 <- function(
  n = 188,
  prob_day = c(64, 64, 60) / 188,
  prob_pot = c(12,12,12,12,12,12, 12,12,12,12, 12,12,12,12,8,12) / 188
){
  # Design 
  day           <- sample(1:3, size = n, replace = TRUE, prob = prob_day)
  pot_treatment <- c(rep(1, 6), rep(2, 4), rep(3, 6))  
  pot_idx       <- sample(1:16, size = n, replace = TRUE, prob = prob_pot)
  treatment     <- pot_treatment[pot_idx]             

  # Parameters 
  kappa_d <- rexp(3, rate = 0.01)        # concentration per day
  tau_t   <- rexp(3, rate = 5)           # pot-level SD per treatment 
  z_p     <- rnorm(16, mean = 0, sd = 1) # standardized per-pot deviate
  eta     <- matrix(rnorm(9, mean = -1, sd = 0.3), nrow = 3, ncol = 3)  

  # per-pot offset
  u_pot <- z_p * tau_t[pot_treatment]    # length 16

  # per-obs 
  eta_i   <- eta[cbind(treatment, day)]  # scalar per obs
  u_i     <- u_pot[pot_idx]              # obs pot offset
  mu_i    <- plogis(eta_i + u_i)         # probability scale
  kappa_i <- kappa_d[day]               
  gammaRC <- rbeta(n, shape1 = mu_i * kappa_i, shape2 = (1 - mu_i) * kappa_i)

  # Return sim data and paramters 
  return(list(
    data = tibble(
      day_idx   = day,
      treatment = treatment,
      pot_idx   = pot_idx,
      gammaRC   = gammaRC
    ),
    parameters = list(
      eta      = eta,        
      kappa_d  = kappa_d,
      tau_t    = tau_t,
      z_p      = z_p,
      u_pot    = u_pot
    )
  ))
}

#### Cell Means Combined Across all estimands ####
dgp_combined_jip <- function(
  n = 188,
  prob_day = c(64, 64, 60) / 188,
  prob_pot = c(12,12,12,12,12,12, 12,12,12,12, 12,12,12,12,8,12) / 188
){
  # Design 
  day           <- sample(1:3, size = n, replace = TRUE, prob = prob_day)
  pot_treatment <- c(rep(1, 6), rep(2, 4), rep(3, 6))  
  pot_idx       <- sample(1:16, size = n, replace = TRUE, prob = prob_pot)
  treatment     <- pot_treatment[pot_idx]

  ## GammaRC Model ## ====================================================
  # Parameters 
  g_kappa_d <- rexp(3, rate = 0.01)        # concentration per day
  g_tau_t   <- rexp(3, rate = 5)           # pot-level SD per treatment 
  g_z_p     <- rnorm(16, mean = 0, sd = 1) # standardized per-pot deviate
  g_eta     <- matrix(rnorm(9, mean = -1, sd = 0.3), nrow = 3, ncol = 3)  

  # per-pot offset
  g_u_pot <- g_z_p * g_tau_t[pot_treatment]    # length 16

  # per-obs 
  g_eta_i   <- g_eta[cbind(treatment, day)]  # scalar per obs
  g_u_i     <- g_u_pot[pot_idx]              # obs pot offset
  g_mu_i    <- plogis(g_eta_i + g_u_i)         # probability scale
  g_kappa_i <- g_kappa_d[day]      
  
  # Draw from Beta Distribution 
  gammaRC <- rbeta(
    n = n,
    shape1 = g_mu_i * g_kappa_i,
    shape2 = (1 - g_mu_i) * g_kappa_i
  )

  ## Phi Model ## ========================================================
  # Parameters
  phi_eta_dt  <- matrix(rnorm(n = 9, mean = 0.8, sd = 0.3),nrow = 3, ncol = 3)
  phi_kappa_d <- rexp(n = 3, rate = 0.001)
  phi_tau_t   <- rexp(n = 3, rate = 3)
  phi_zp      <- rnorm(n = 16, mean = 0, sd = 1)

  # Per pot offect 
  phi_up <- phi_tau_t[pot_treatment] * phi_zp

  # Per observation 
  phi_eta_dti   <- phi_eta_dt[cbind(treatment, day)] 
  phi_kappa_d_i <- phi_kappa_d[day] 
  phi_up_i      <- phi_up[pot_idx] 
  phi_mu_i      <- plogis(phi_eta_dti + phi_up_i)

  # Draw from beta distribution 
  phi_po <- rbeta(
    n = n,
    shape1 = phi_mu_i * phi_kappa_d_i,
    shape2 = (1 - phi_mu_i) * phi_kappa_d_i
  )

  ## Psi Model ## ===========================================================
  # Paramters 
  psi_eta_dt  <- matrix(rnorm(n = 9, mean = 0.7, sd = 0.3), nrow = 3, ncol = 3)
  psi_kappa_d <- rexp(n = 3, rate = 0.01) 
  psi_tau_t   <- rexp(n = 3, rate = 3)
  psi_zp      <- rnorm(n = 16, mean = 0, sd = 1)

  # Per pot offcet 
  psi_up <- psi_tau_t[pot_treatment] * psi_zp

  # Per obs 
  psi_eta_td_i  <- psi_eta_dt[cbind(treatment, day)]
  psi_up_i      <- psi_up[pot_idx]
  psi_kappa_d_i <- psi_kappa_d[day] 
  psi_mu_i      <- plogis(psi_eta_td_i + psi_up_i)

  # Draw from beta Distribution 
  psi_Eo <- rbeta(
    n = n,
    shape1 = psi_mu_i * psi_kappa_d_i,
    shape2 = (1 - psi_mu_i) * psi_kappa_d_i
  )

  ## Ro Model ## ============================================================
  # Paramters 
  ro_eta_dt   <- matrix(rnorm(n = 9, mean = 0.1, sd = 0.7),nrow = 3, ncol = 3)
  ro_kappa_d  <- rexp(n = 3, rate = 0.01) 
  ro_tau_t    <- rexp(n = 3, rate = 2)
  ro_zp       <- rnorm(n = 16, mean = 0, sd = 1)
  
  # Pot Offcet 
  ro_up <- ro_tau_t[pot_treatment] * ro_zp

  # Per obs
  ro_eta_dt_i  <- ro_eta_dt[cbind(treatment, day)]
  ro_kappa_d_i <- ro_kappa_d[day]
  ro_up_i      <- ro_up[pot_idx]
  ro_mu_i      <- plogis(ro_eta_dt_i + ro_up_i)

  # Draw from Beta Distribution 
  delta_Ro <- rbeta(
    n = n,
    shape1 = ro_mu_i * ro_kappa_d_i,
    shape2 = (1 - ro_mu_i) * ro_kappa_d_i
  )

  ## Combine in one dataset ##
  sim_data <- data.frame(
    day_idx   = day,
    treatment = treatment,
    pot_idx   = pot_idx,
    delta_Ro = delta_Ro,
    psi_Eo   = psi_Eo,
    phi_po   = phi_po,
    gammaRC  = gammaRC
  )

  # Combine the paramters in a list 
  paramters <- list(
    gammaRC = list(
      kappa_d = g_kappa_d,   
      tau     = g_tau_t,          
      zp      = g_z_p,
      eta     = g_eta  
    ),
    phi_Po = list(
      eta   = phi_eta_dt,
      kappa = phi_kappa_d,
      tau   = phi_tau_t,
      zp    = phi_zp
    ),
    psi_Eo = list(
      eta   = psi_eta_dt,
      kappa = psi_kappa_d,
      tau   = psi_tau_t,
      zp    = psi_zp
    ),
    delta_Ro = list(
      eta   = ro_eta_dt,
      kappa = ro_kappa_d,
      tau   = ro_tau_t,
      zp    = ro_zp
    )
  )
  # Return Simulated Data and Parameters 
  return(list(data = sim_data, parameters = paramters))
}
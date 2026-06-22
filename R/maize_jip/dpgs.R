
#### Function to for Simulating Data Generative Processes to test the model ####

#### Cell- Means model Without Pooling ####
dgp_gammaRC_jip_A <- function(
  n = 188,
  prob_treatment = c(68, 72, 48) / 188,
  prob_day       = c(64, 64, 60) / 188
){
  # Sumulate Day and  Treatment
  day       <- sample(1:3, size = n, replace = TRUE, prob = prob_day)
  treatment <- sample(1:3, size = n, replace = TRUE, prob = prob_treatment)

  # Simulate Prior Matrix from Normal Distribution for each d and t 
  eta <- matrix(rnorm(n = 9, mean = -1, sd = 1),nrow = 3,ncol = 3)
  mu_mat <- plogis(eta)
  
  # For every day D allow seperate kappa variation 
  kappa_d <- rexp(n = 3, rate = 0.01)

  # Per obervation pair
  mu_i    <- mu_mat[treatment, day]
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
      day = day,
      treatment = treatment,
      gamma_RC = gammaRC
    ),
    parameters = tibble(
      day = 1:3,
      kappa_per_day = kappa_d
    ),
    cell_means_logit = eta, 
    cell_means_prob  = mu_mat
  ))
}

#### Cell-Means Model With Pot Pooling ####
dgp_gammaRC_jip_B <- function(
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
  tau_t   <- rexp(3, rate = 1)           # pot-level SD per treatment 
  z_p     <- rnorm(16, mean = 0, sd = 1) # standardized per-pot deviate
  eta     <- matrix(rnorm(9, mean = -1, sd = 1), nrow = 3, ncol = 3)  

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
      day       = day,
      treatment = treatment,
      pot       = pot_idx,
      gamma_RC  = gammaRC
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
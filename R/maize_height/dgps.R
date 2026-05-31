
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

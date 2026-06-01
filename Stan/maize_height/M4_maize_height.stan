
// Data Input Block
data{
    int<lower=1> N;
    int<lower=1> G;                           // Grid Size for ploting  
    int<lower=0,upper=1> prior_only;
    array[N]int<lower=1,upper=3> treatment;
    array[N]int<lower=1,upper=4> day_idx;     // Day index
    vector[N] t;                              // Actual time 
    vector[N] plant_height;
}
// Model Paramters 
parameters {
    vector[3] log_A;
    vector[3] log_K;
    vector[3] tau_star;
    vector<lower=0>[4] sigma;
}
// Parameter Transformation 
transformed parameters {
    vector[3] A = exp(log_A);
    vector[3] K = exp(log_K);
}

// Model 
model{
    // Priors
    log_A    ~ normal(log(40), 0.3);
    log_K    ~ normal(log(0.15), 0.5);
    tau_star ~ normal(9, 3);
    sigma    ~ exponential(2);

    // Model Likelihood 
    if (prior_only == 0) {
        for (i in 1:N) {
            real mu = A[treatment[i]] * exp(-exp(-K[treatment[i]] * (t[i] - tau_star[treatment[i]])));
            plant_height[i] ~ lognormal(log(mu), sigma[day_idx[i]]);
        }
    }
}

// Additional Calculations 
generated quantities {
    
    // Pairwise contrasts 
    real log_A_emf_vs_control  = log_A[3] - log_A[1];
    real log_A_emf_vs_sham     = log_A[3] - log_A[2];
    real log_A_sham_vs_control = log_A[2] - log_A[1];

    real log_K_emf_vs_control  = log_K[3] - log_K[1];
    real log_K_emf_vs_sham     = log_K[3] - log_K[2];
    real log_K_sham_vs_control = log_K[2] - log_K[1];

    // tau_star contrasts on the raw scale units: days
    real tau_emf_vs_control  = tau_star[3] - tau_star[1];
    real tau_emf_vs_sham     = tau_star[3] - tau_star[2];
    real tau_sham_vs_control = tau_star[2] - tau_star[1];

    // Fitted median curves on a fine grid for plotting 
    // 51 points from day 0 to day 35
    vector[G] t_grid;
    matrix[3, G] mu_grid;
    for (g in 1:G) {
        t_grid[g] = (g - 1) * 35.0 / (G - 1);   
        for (tr in 1:3) {
            mu_grid[tr, g] = A[tr] * exp(-exp(-K[tr] * (t_grid[g] - tau_star[tr])));
        }
    }
    // for loo and PPC 
    vector[N] log_lik;
    array[N] real h_rep;
    for (i in 1:N) {
        real mu_i = A[treatment[i]] * exp(-exp(-K[treatment[i]] * (t[i] - tau_star[treatment[i]])));
        log_lik[i] = lognormal_lpdf(plant_height[i] | log(mu_i), sigma[day_idx[i]]);
        h_rep[i]   = lognormal_rng(log(mu_i), sigma[day_idx[i]]);
    }

}

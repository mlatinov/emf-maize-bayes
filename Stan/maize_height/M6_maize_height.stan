// Data Input Block 
data{
    // Obervations in the data N and Prior only Switch
    int<lower=1> N;
    int<lower=0,upper=1> prior_only;
    
    // Day Variables 
    array[N] int<lower=1,upper=4> day_idx; // Day Index
    vector[N] t;                        // Actual time 

    // Pot Variables               
    array[N]  int<lower=1,upper=16> pot_idx; // Pot Index
    array[16] int<lower=1, upper=3> pot_treatment;

    // Target Plant height
    vector[N] plant_height;
} 
// Model Paramters Block 
parameters{
    // Population Level Paramters per treatment 
    vector[3] a_t; // asymptote 
    vector[3] b_t; // rate level 
    vector[3] c_t; // inflection point 

    // Pot-scatter SD per treatment
    vector<lower=0>[3] xi_A_t;
    vector<lower=0>[3] xi_K_t;
    vector<lower=0>[3] xi_tau_t;

    // Unit-scale pot deviation
    vector[16] z_A;
    vector[16] z_K;
    vector[16] z_tau;

    // Residual SD per day 
    vector<lower=0.001>[4] sigma_d;
}
// Paramters Transformations 
transformed parameters{
    // Non centered Paramterarization 
    vector[16] A_p;
    vector[16] K_p;
    vector[16] tau_p;
    for(i in 1:16){
        A_p[i] = exp(a_t[pot_treatment[i]] + xi_A_t[pot_treatment[i]] * z_A[i]);
        K_p[i] = exp(b_t[pot_treatment[i]] + xi_K_t[pot_treatment[i]] * z_K[i]);
        tau_p[i] = c_t[pot_treatment[i]] + xi_tau_t[pot_treatment[i]] * z_tau[i];
    }
}
// Model 
model{
    // Priors //

    // Population Level 
    a_t ~ normal(log(40), 0.3);
    b_t ~ normal(log(0.15), 0.5);
    c_t ~ normal(9, 3);

    // Pot-scatter SD 
    xi_A_t ~ exponential(10);
    xi_K_t ~ exponential(10);
    xi_tau_t ~ exponential(10);

    // Unit-scale pot deviation
    z_A ~ normal(0, 1);
    z_K ~ normal(0, 1);
    z_tau ~ normal(0, 1);

    // Residual SD
    sigma_d ~ exponential(2);

    // Model Likelihood //
    if(prior_only == 0){
        for(i in 1:N){
            // Golmez Equantion to Calculate mu
            real mu = A_p[pot_idx[i]] * exp(-exp(-K_p[pot_idx[i]] * (t[i] - tau_p[pot_idx[i]])));
            // Plant Height with Sigma per day index 
            plant_height[i] ~ lognormal(log(mu), sigma_d[day_idx[i]]);
        }
    }
}
// Additional Calculations 
generated quantities {
    
    // Population-level natural-scale parameters 
    vector[3] A_pop = exp(a_t);
    vector[3] K_pop = exp(b_t);
    
    // Contrasts on population curve parameters
    real log_A_emf_vs_control  = a_t[3] - a_t[1];
    real log_A_emf_vs_sham     = a_t[3] - a_t[2];
    real log_A_sham_vs_control = a_t[2] - a_t[1];
    
    // log-rate contrasts
    real log_K_emf_vs_control  = b_t[3] - b_t[1];
    real log_K_emf_vs_sham     = b_t[3] - b_t[2];
    real log_K_sham_vs_control = b_t[2] - b_t[1];
    
    // Inflection contrasts 
    real tau_emf_vs_control  = c_t[3] - c_t[1];
    real tau_emf_vs_sham     = c_t[3] - c_t[2];
    real tau_sham_vs_control = c_t[2] - c_t[1];
    
    // For loo and PPC
    vector[N] log_lik;
    array[N] real h_rep;
    for (i in 1:N) {
        real mu = A_p[pot_idx[i]] * exp(-exp(-K_p[pot_idx[i]] * (t[i] - tau_p[pot_idx[i]])));
        log_lik[i] = lognormal_lpdf(plant_height[i] | log(mu), sigma_d[day_idx[i]]);
        h_rep[i]   = lognormal_rng(log(mu), sigma_d[day_idx[i]]);
    }

}

// Function Block 
functions{
    #include ../lib/effect_probabilities.stanfunctions
    #include ../lib/diagnostics.stanfunctions
}

// Data Input Block 
data{
    // Obervations in the data N and Prior only Switch
    int<lower=1> N;
    int<lower=0,upper=1> prior_only;
    int<lower=1> G; // Grid points for ploting 
    
    // Day Variables 
    array[N] int<lower=1,upper=4> day_idx; // Day Index
    vector[N] t;                           // Actual time 

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
    c_t ~ normal(9, 1);

    // Pot-scatter SD 
    xi_A_t ~ exponential(10);
    xi_K_t ~ exponential(10);
    xi_tau_t ~ exponential(3);

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

    // === Population-level natural-scale parameters ========================
    vector[3] A_pop = exp(a_t);
    vector[3] K_pop = exp(b_t);
    
    // === Shared linear predictors ================================================
    vector[N] mu_full;          // pot-level prediction (cm scale)
    vector[N] mu_pop;           // population-level prediction (cm scale)  
    vector[N] log_mu_full;      // log-scale, for residual diagnostics & loo
    vector[N] log_mu_pop;       // log-scale equivalent
    vector[N] sigma_per_obs;
    for (i in 1:N) {
        // Pot-level curve
        mu_full[i] = A_p[pot_idx[i]] * exp(-exp(-K_p[pot_idx[i]] * (t[i] - tau_p[pot_idx[i]])));
        // Population-level curve (treatment lookup via pot_treatment)
        {
            int tr = pot_treatment[pot_idx[i]];
            mu_pop[i] = A_pop[tr] * exp(-exp(-K_pop[tr] * (t[i] - c_t[tr])));
        }
        log_mu_full[i]   = log(mu_full[i]);
        log_mu_pop[i]    = log(mu_pop[i]);
        sigma_per_obs[i] = sigma_d[day_idx[i]];
    }

    // === Pairwise contrasts ===============================================
    real log_A_emf_vs_control  = a_t[3] - a_t[1];
    real log_A_emf_vs_sham     = a_t[3] - a_t[2];
    real log_A_sham_vs_control = a_t[2] - a_t[1];
    
    real log_K_emf_vs_control  = b_t[3] - b_t[1];
    real log_K_emf_vs_sham     = b_t[3] - b_t[2];
    real log_K_sham_vs_control = b_t[2] - b_t[1];
    
    real tau_emf_vs_control  = c_t[3] - c_t[1];
    real tau_emf_vs_sham     = c_t[3] - c_t[2];
    real tau_sham_vs_control = c_t[2] - c_t[1];

    // === Probability-of-direction indicators ==============================
    int prob_log_A_emf_vs_control_positive  = is_positive(log_A_emf_vs_control);
    int prob_log_A_emf_vs_sham_positive     = is_positive(log_A_emf_vs_sham);
    int prob_log_A_sham_vs_control_positive = is_positive(log_A_sham_vs_control);

    int prob_log_K_emf_vs_control_positive  = is_positive(log_K_emf_vs_control);
    int prob_log_K_emf_vs_sham_positive     = is_positive(log_K_emf_vs_sham);
    int prob_log_K_sham_vs_control_positive = is_positive(log_K_sham_vs_control);
    
    int prob_tau_emf_vs_control_positive    = is_positive(tau_emf_vs_control);
    int prob_tau_emf_vs_sham_positive       = is_positive(tau_emf_vs_sham);
    int prob_tau_sham_vs_control_positive   = is_positive(tau_sham_vs_control);

    // === ROPE indicators ==================================================
    int rope_log_A_emf_vs_control  = in_rope(log_A_emf_vs_control,  -0.05, 0.05);
    int rope_log_A_emf_vs_sham     = in_rope(log_A_emf_vs_sham,     -0.05, 0.05);
    int rope_log_A_sham_vs_control = in_rope(log_A_sham_vs_control, -0.05, 0.05);
    
    int rope_log_K_emf_vs_control  = in_rope(log_K_emf_vs_control,  -0.10, 0.10);
    int rope_log_K_emf_vs_sham     = in_rope(log_K_emf_vs_sham,     -0.10, 0.10);
    int rope_log_K_sham_vs_control = in_rope(log_K_sham_vs_control, -0.10, 0.10);
    
    int rope_tau_emf_vs_control    = in_rope(tau_emf_vs_control,    -0.5,  0.5);
    int rope_tau_emf_vs_sham       = in_rope(tau_emf_vs_sham,       -0.5,  0.5);
    int rope_tau_sham_vs_control   = in_rope(tau_sham_vs_control,   -0.5,  0.5);

    // === Variance decomposition ======================================
    real var_pot_contribution = variance(log_mu_full - log_mu_pop);
    real var_res_full = mean(square(sigma_d));
    real var_res_pop  = mean(square(sigma_d)) + var_pot_contribution;

    real bayes_R2_full       = bayes_R2_general(log_mu_full, var_res_full);
    real bayes_R2_population = bayes_R2_general(log_mu_pop,  var_res_pop);
    real bayes_R2_pot_contribution = bayes_R2_full - bayes_R2_population;
    
    // Residual diagnostics ===============================================================
    vector[N] log_resid_pearson;
    vector[N] pit;
    if (prior_only == 0) {
        log_resid_pearson = pearson_residuals_lognorm(plant_height, mu_full, sigma_per_obs);
        pit               = lognormal_pit(plant_height, mu_full, sigma_per_obs);
    } else {
        log_resid_pearson = rep_vector(0.0, N);
        pit               = rep_vector(0.5, N);  
    }
    
    // === Curve prediction grid ======================
    vector[G] t_grid;
    matrix[3, G] mu_pop_grid;
    for (g in 1:G) {
        t_grid[g] = (g - 1) * 35.0 / (G - 1);
        for (tr in 1:3) {
            mu_pop_grid[tr, g] = A_pop[tr] * exp(-exp(-K_pop[tr] * (t_grid[g] - c_t[tr])));
        }
    }

    // === Population median at each measurement day ========================
    matrix[3, 4] median_at_day;
    {
        array[4] real day_actual = {7, 14, 20, 27};
        for (tr in 1:3) for (d in 1:4) {
            median_at_day[tr, d] = A_pop[tr] * 
                                    exp(-exp(-K_pop[tr] * (day_actual[d] - c_t[tr])));
        }
    }

    // === Likelihood for LOO and replicates for PPC ========================
    vector[N] log_lik;
    array[N] real h_rep;
    for (i in 1:N) {
        log_lik[i] = lognormal_lpdf(plant_height[i] | log_mu_full[i], sigma_per_obs[i]);
        h_rep[i]   = lognormal_rng(log_mu_full[i], sigma_per_obs[i]);
    }
}

// Function Block 
functions{
    #include ../lib/effect_probabilities.stanfunctions
}
// Data input Block 
data{
    int<lower=1> N; 
    int<lower=0,upper=1> prior_only;
    vector[N] plant_height;
    array[16]int<lower=1, upper=3> pot_treatment;
    array[N] int<lower=1, upper=3> treatment;
    array[N] int<lower=1, upper=4> day_idx;
    array[N] int<lower=1, upper=16> pot_idx;
}

// Model paramters 
parameters{
    matrix[3, 4] alpha_dt;            // log-median height of treatment on day 
    vector[16]   z_bar_p;             // Unit-scale pot deviation
    vector<lower=0.001>[3] zeta_t;    // SD of pot offsets within treatment 
    vector<lower=0.001>[4] sigma_d;   // residual SD on day 
}

// Transform Parameters Block 
transformed parameters{
    // pot-level offsets (one per pot, p=1,…,16)
    vector[16] zp;               
    for (p in 1:16) {
        zp[p] = zeta_t[pot_treatment[p]] * z_bar_p[p];
    }
}

// Model Block 
model{
    // Priors 
    to_vector(alpha_dt) ~ normal(3, 1);
    z_bar_p ~ normal(0, 1);
    zeta_t  ~ exponential(10);
    sigma_d    ~ exponential(2); 

    // Model Likelihood 
    if(prior_only == 0){
        for(i in 1:N){
            real mu_i = alpha_dt[treatment[i], day_idx[i]] + zp[pot_idx[i]];
            plant_height[i] ~ lognormal(mu_i, sigma_d[day_idx[i]]);
        }
    }
}

// Additional Calculations 
generated quantities {

    // Cell-level summaries ==================================================
    // Per-day cell medians on the cm scale (for plotting / interpretation)
    matrix[3, 4] median_height;
    for (t in 1:3) for (d in 1:4) median_height[t, d] = exp(alpha_dt[t, d]);

    // Pairwise contrasts per day =====================================================
    vector[4] log_diff_emf_vs_control;
    vector[4] log_diff_emf_vs_sham;
    vector[4] log_diff_sham_vs_control;
    for (d in 1:4) {
        log_diff_emf_vs_control[d]  = alpha_dt[3, d] - alpha_dt[1, d];
        log_diff_emf_vs_sham[d]     = alpha_dt[3, d] - alpha_dt[2, d];
        log_diff_sham_vs_control[d] = alpha_dt[2, d] - alpha_dt[1, d];
    }

    // Probability-of-direction indicators pre Day =====================================
    vector[4] prob_emf_control_positive;
    vector[4] prob_emf_sham_positive;
    vector[4] prob_sham_control_positive;
    for(d in 1:4) {
        prob_emf_control_positive[d] = is_positive(log_diff_emf_vs_control[d]);
        prob_emf_sham_positive[d]    = is_positive(log_diff_emf_vs_sham[d]);
        log_diff_sham_vs_control[d]  = is_positive(log_diff_sham_vs_control[d]);
    }

    // ROPE indicators ===============================================================
    vector[4] rope_emf_control_diff;
    vector[4] rope_emf_sham_diff;
    vector[4] rope_sham_control_diff;
    for(d in 1:4) {
        rope_emf_control_diff[d]  = in_rope(log_diff_emf_vs_control[d], -0.05, 0.05);
        rope_emf_sham_diff[d]     = in_rope(log_diff_emf_vs_sham[d],    -0.05, 0.05);
        rope_sham_control_diff[d] = in_rope(log_diff_sham_vs_control[d] -0.05, 0.05);
    }

    // Variance decomposition ======================================================

    // For loo and PPC
    vector[N] log_lik;
    array[N] real h_rep;
    for (i in 1:N) {
        real mu_i = alpha_dt[treatment[i], day_idx[i]] + zp[pot_idx[i]];
        log_lik[i] = lognormal_lpdf(plant_height[i] | mu_i, sigma_d[day_idx[i]]);
        h_rep[i]   = lognormal_rng(mu_i, sigma_d[day_idx[i]]);
    }

}

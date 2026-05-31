
// Input Data block 
data{
    int<lower=1> N;
    vector[N] plant_height;
    array[N] int<lower=1,upper=16> pot_id;
    array[N] int<lower=1,upper=4> day;
    int<lower=0,upper=1> prior_only;
}

// Model Paramters 
parameters{
    matrix[16, 4] alpha;
    vector<lower=0.001>[4] sigma;
}

// Model 
model{
    // Priors 
    to_vector(alpha) ~ normal(3, 1);
    sigma            ~ exponential(2);

    // Model likelihood 
    if(!prior_only){
        for(i in 1:N){
            plant_height[i] ~ lognormal(alpha[pot_id[i], day[i]], sigma[day[i]]);
        }
    }
}

// Additional Calculations 
generated quantities {
    // Per-pot-per-day medians on cm scale
    matrix[16, 4] median_height;
    for (p in 1:16) for (d in 1:4) median_height[p, d] = exp(alpha[p, d]);

    // Recover treatment-level means by averaging the relevant pot rows
    matrix[3, 4] treatment_mean_log;
    for (d in 1:4) {
        treatment_mean_log[1, d] = mean(alpha[1:6,   d]);   // Control: pots 1-6
        treatment_mean_log[2, d] = mean(alpha[7:10,  d]);   // Sham:    pots 7-10
        treatment_mean_log[3, d] = mean(alpha[11:16, d]);   // EMF:     pots 11-16
    }

    // Per-day treatment contrasts (log scale)
    vector[4] log_diff_emf_vs_control;
    vector[4] log_diff_emf_vs_sham;
    vector[4] log_diff_sham_vs_control;
    for (d in 1:4) {
        log_diff_emf_vs_control[d]  = treatment_mean_log[3, d] - treatment_mean_log[1, d];
        log_diff_emf_vs_sham[d]     = treatment_mean_log[3, d] - treatment_mean_log[2, d];
        log_diff_sham_vs_control[d] = treatment_mean_log[2, d] - treatment_mean_log[1, d];
    }

    // For loo and PPC
    vector[N] log_lik;
    array[N] real h_rep;
    for (i in 1:N) {
        log_lik[i] = lognormal_lpdf(plant_height[i] | alpha[pot_id[i], day[i]], sigma[day[i]]);
        h_rep[i]   = lognormal_rng(alpha[pot_id[i], day[i]], sigma[day[i]]);
    }
}

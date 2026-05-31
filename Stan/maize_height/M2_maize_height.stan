
// Data input Block 
data{
    int<lower=1> N;
    array[N] int<lower=1, upper=3> treatment;
    array[N] int<lower=1, upper=4> day;
    vector[N] plant_height;
    int<lower=0,upper=1> prior_only;
}
// Model Paramters 
parameters{
    matrix[3, 4] alpha;
    vector<lower=0.001>[4] sigma;
}
// Model 
model{
    // Priors 
    to_vector(alpha) ~ normal(3, 1);
    sigma            ~ exponential(2);

    // Model Likelihood
    if(prior_only == 0){
        for(i in 1:N){
            plant_height[i] ~ lognormal(alpha[treatment[i],day[i]], sigma[day[i]]); 
        }
    }
}
// Additional Calculations 
generated quantities {
    // Group-level medians on the cm scale 
    matrix[3, 4] median_height;
    for (t in 1:3) for (d in 1:4) median_height[t, d] = exp(alpha[t, d]);

    // Per-day log-scale contrasts treatment effects at each day
    vector[4] log_diff_emf_vs_control;
    vector[4] log_diff_emf_vs_sham;
    vector[4] log_diff_sham_vs_control;
    for (d in 1:4) {
        log_diff_emf_vs_control[d]  = alpha[3, d] - alpha[1, d];
        log_diff_emf_vs_sham[d]     = alpha[3, d] - alpha[2, d];
        log_diff_sham_vs_control[d] = alpha[2, d] - alpha[1, d];
    }

    // Per-day ratio contrasts 
    vector[4] ratio_emf_vs_control = exp(log_diff_emf_vs_control);
    vector[4] ratio_emf_vs_sham    = exp(log_diff_emf_vs_sham);

    // For loo and PPC
    vector[N] log_lik;
    array[N] real h_rep;
    for (i in 1:N) {
        log_lik[i] = lognormal_lpdf(plant_height[i] | alpha[treatment[i], day[i]], sigma[day[i]]);
        h_rep[i]   = lognormal_rng(alpha[treatment[i], day[i]], sigma[day[i]]);
    }
}

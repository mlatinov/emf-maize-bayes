// Function Block 
functions{
    #include ../lib/effect_probabilities.stanfunctions
    #include ../lib/diagnostics.stanfunctions
}
// Data Input Block
data{
    int<lower = 1> N;
    int<lower = 0, upper = 1> prior_only;
    vector<lower = 0, upper =1>[N] gammaRC;
    array[N] int<lower=1, upper=3> day_idx;
    array[N] int<lower=1, upper=3> treatment;
}
// Model Paramaters 
parameters{
    matrix[3,3] eta;
    vector<lower = 0>[3] kappa_d;
}
// Transform Model Paramters 
transformed parameters {
   matrix[3, 3] mu = inv_logit(eta);
}
// Model Block 
model{
    // Priors 
    to_vector(eta) ~ normal(-1, 1);
    kappa_d ~ exponential(0.01);

    // Model Likelihood 
    if(prior_only == 0){
        for(i in 1:N){
            gammaRC[i] ~ beta(
                mu[treatment[i], day_idx[i]] * kappa_d[day_idx[i]],
                (1 - mu[treatment[i], day_idx[i]]) * kappa_d[day_idx[i]]
            );
        }
    }
}
// Additional Calculations 
generated quantities {
    // Posterior or prior predictive replicates
    vector[N] gammaRC_rep;
    for (i in 1:N) {
        real m = mu[treatment[i], day_idx[i]];
        real k = kappa_d[day_idx[i]];
        gammaRC_rep[i] = beta_rng(m * k, (1 - m) * k);
    }
    // Derived RC/ABS per cell 
    matrix[3,3] RC_ABS;
    for (t in 1:3)
        for (d in 1:3)
            RC_ABS[t,d] = mu[t,d] / (1 - mu[t,d]);

}


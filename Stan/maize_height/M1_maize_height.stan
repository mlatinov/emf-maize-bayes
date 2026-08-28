
// Data Input Block 
data{
    int <lower = 1> N; // Number of observations
    array[N]int<lower = 1, upper = 3>  treatment;
    vector<lower = 0>[N] plant_height;
    int<lower = 0, upper = 1>  prior_only;
}
// Model Parameters block
parameters{
    vector[3] alpha;
    real<lower=0.001> sigma;
}
// Model 
model{
    // Priors 
    alpha ~ normal(3, 1);
    sigma ~ exponential(2);

    // Model Likelihood 
    if(prior_only == 0){
        for(i in 1:N){
            plant_height[i] ~ lognormal(alpha[treatment[i]], sigma);
        }
    }
}
// Additional Calculations 
generated quantities {

  // Group-level posteriors 
  real median_control = exp(alpha[1]);
  real median_sham    = exp(alpha[2]);
  real median_emf     = exp(alpha[3]);

  // Log-scale contrasts 
  real log_diff_emf_vs_control  = alpha[3] - alpha[1];
  real log_diff_emf_vs_sham     = alpha[3] - alpha[2];
  real log_diff_sham_vs_control = alpha[2] - alpha[1];

  // Ratio-scale contrasts 
  real ratio_emf_vs_control = exp(alpha[3] - alpha[1]);
  real ratio_emf_vs_sham    = exp(alpha[3] - alpha[2]);

  // PPC
  vector[N] log_lik;
  array[N] real h_rep;
  for (i in 1:N) {
    log_lik[i] = lognormal_lpdf(plant_height[i] | alpha[treatment[i]], sigma);
    h_rep[i]   = lognormal_rng(alpha[treatment[i]], sigma);
  }

}

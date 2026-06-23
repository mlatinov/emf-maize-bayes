// Function Block 
functions{
    #include ../lib/effect_probabilities.stanfunctions
    #include ../lib/diagnostics.stanfunctions
}
// Data input Block 
data{
    int<lower = 1> N;
    int<lower = 0> prior_only;
    vector<lower = 0, upper = 1>[N] gammaRC; 
    array[N] int<lower = 1, upper = 16> pot_idx;
    array[N] int<lower = 1, upper = 3> treatment;
    array[N] int<lower = 1, upper = 3> day_idx;
    array[16]int<lower = 1, upper = 3> pot_treatment;
}
// Model Paramters 
parameters{
    matrix[3, 3] eta;              // Cell-Means Matrix Day x Treatment
    vector<lower = 0>[3] kappa_d;  // Range of Dispersion  
    vector<lower = 0>[3] tau_t;    // Within Treatment Dispersion
    vector[16] zp;                 // Scaling Parameter per pot 
}
// Paramter Transformations 
transformed parameters {
    // Calculate per pot dispersion 
    vector[16] u_p;
    for (p in 1:16)
        u_p[p] = zp[p] * tau_t[pot_treatment[p]];
        
    // Calculate per observation mean on the prob scale  
    vector[N] mu;
    for (i in 1:N)
        mu[i] = inv_logit(eta[treatment[i], day_idx[i]] + u_p[pot_idx[i]]);
}
// Model 
model{
    // Priors
    to_vector(eta) ~ normal(-1, 1);
    kappa_d        ~ exponential(0.01);
    tau_t          ~ exponential(2);
    zp             ~ normal(0,1);
    
    // Model Likelihood
    if (prior_only == 0){
        for (i in 1:N) {
            real k = kappa_d[day_idx[i]];
            gammaRC[i] ~ beta( mu[i] * k, (1 - mu[i]) * k );
        }
    }
}
// Additional Calculations 
generated quantities {

    // === Population-level (pot-marginal) cell means =======================
    matrix[3,3] mu_pop;
    for (t in 1:3)
        for (d in 1:3)
            mu_pop[t,d] = inv_logit(eta[t,d]);   // probability scale

    // Derived odds quantity RC/ABS = mu / (1 - mu), per cell
    matrix[3,3] RC_ABS;
    for (t in 1:3)
        for (d in 1:3)
            RC_ABS[t,d] = mu_pop[t,d] / (1 - mu_pop[t,d]);

    // === Pairwise contrasts, per day =====================================
    matrix[3,3] logit_emf_vs_control;   // day reused as vector; store per day
    vector[3] dlogit_emf_vs_control;
    vector[3] dlogit_emf_vs_sham;
    vector[3] dlogit_sham_vs_control;
    vector[3] dprob_emf_vs_control;
    vector[3] dprob_emf_vs_sham;
    vector[3] dprob_sham_vs_control;
    for (d in 1:3) {
        dlogit_emf_vs_control[d]  = eta[3,d] - eta[1,d];
        dlogit_emf_vs_sham[d]     = eta[3,d] - eta[2,d];
        dlogit_sham_vs_control[d] = eta[2,d] - eta[1,d];

        dprob_emf_vs_control[d]   = mu_pop[3,d] - mu_pop[1,d];
        dprob_emf_vs_sham[d]      = mu_pop[3,d] - mu_pop[2,d];
        dprob_sham_vs_control[d]  = mu_pop[2,d] - mu_pop[1,d];
    }

    // === Probability-of-direction indicators (per day) ====================
    array[3] int pd_emf_vs_control;
    array[3] int pd_emf_vs_sham;
    array[3] int pd_sham_vs_control;
    for (d in 1:3) {
        pd_emf_vs_control[d]  = is_positive(dlogit_emf_vs_control[d]);
        pd_emf_vs_sham[d]     = is_positive(dlogit_emf_vs_sham[d]);
        pd_sham_vs_control[d] = is_positive(dlogit_sham_vs_control[d]);
    }

    // === ROPE indicators ==================================================
    real rope_prob_half = 0.02;
    real rope_logit_half = 0.10;
    array[3] int rope_prob_emf_vs_control;
    array[3] int rope_prob_emf_vs_sham;
    array[3] int rope_prob_sham_vs_control;
    array[3] int rope_logit_emf_vs_control;
    array[3] int rope_logit_emf_vs_sham;
    array[3] int rope_logit_sham_vs_control;
    for (d in 1:3) {
        rope_prob_emf_vs_control[d]  = in_rope(dprob_emf_vs_control[d],  -rope_prob_half, rope_prob_half);
        rope_prob_emf_vs_sham[d]     = in_rope(dprob_emf_vs_sham[d],     -rope_prob_half, rope_prob_half);
        rope_prob_sham_vs_control[d] = in_rope(dprob_sham_vs_control[d], -rope_prob_half, rope_prob_half);
        rope_logit_emf_vs_control[d]  = in_rope(dlogit_emf_vs_control[d],  -rope_logit_half, rope_logit_half);
        rope_logit_emf_vs_sham[d]     = in_rope(dlogit_emf_vs_sham[d],     -rope_logit_half, rope_logit_half);
        rope_logit_sham_vs_control[d] = in_rope(dlogit_sham_vs_control[d], -rope_logit_half, rope_logit_half);
    }

    // === Pot-level heterogeneity summaries ================================
    real sd_u_p = sd(u_p); // overall realized pot-offset SD

    // === Variance decomposition / Bayes R2 ================================
    vector[N] eta_full;      // cell mean + pot offset 
    vector[N] eta_pop;       // cell mean only 
    for (i in 1:N) {
        eta_full[i] = eta[treatment[i], day_idx[i]] + u_p[pot_idx[i]];
        eta_pop[i]  = eta[treatment[i], day_idx[i]];
    }
    // Residual variance on the logit scale: Beta variance mapped approx via
    real var_pot_contribution = variance(eta_full - eta_pop);

    // Mean Beta variance across observations, pushed to logit scale via
    // first-order delta method: Var(logit) ~ Var(mu) / (mu(1-mu))^2.
    vector[N] beta_var_obs;
    for (i in 1:N) {
        real m = inv_logit(eta_full[i]);
        real k = kappa_d[day_idx[i]];
        real v = m * (1 - m) / (k + 1);             // Beta variance on prob scale
        beta_var_obs[i] = v / square(m * (1 - m)); //  logit scale (delta method)
    }
    real var_res_full = mean(beta_var_obs);
    real var_res_pop  = mean(beta_var_obs) + var_pot_contribution;

    real bayes_R2_full       = bayes_R2_general(eta_full, var_res_full);
    real bayes_R2_population = bayes_R2_general(eta_pop,  var_res_pop);
    real bayes_R2_pot_contribution = bayes_R2_full - bayes_R2_population;

    // === Residual diagnostics =============================================
    vector[N] resid_pearson;
    vector[N] pit;
    if (prior_only == 0) {
        for (i in 1:N) {
            real m = inv_logit(eta_full[i]);
            real k = kappa_d[day_idx[i]];
            real a = m * k;
            real b = (1 - m) * k;
            real sd_i = sqrt(m * (1 - m) / (k + 1));
            resid_pearson[i] = (gammaRC[i] - m) / sd_i;
            pit[i]           = beta_cdf(gammaRC[i] | a, b);
        }
    } else {
        resid_pearson = rep_vector(0.0, N);
        pit           = rep_vector(0.5, N);
    }

    // === Likelihood for LOO and replicates for PPC ========================
    vector[N] log_lik;
    array[N] real gammaRC_rep;
    for (i in 1:N) {
        real m = inv_logit(eta_full[i]);
        real k = kappa_d[day_idx[i]];
        real a = m * k;
        real b = (1 - m) * k;
        log_lik[i]     = beta_lpdf(gammaRC[i] | a, b);
        gammaRC_rep[i] = beta_rng(a, b);
    }
}

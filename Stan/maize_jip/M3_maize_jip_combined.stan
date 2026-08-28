// Function Block
functions{
    #include ../lib/effect_probabilities.stanfunctions
    #include ../lib/diagnostics.stanfunctions
}
// Input Data
data{
    // Settings ====================================
    int<lower = 0> N;
    int<lower = 0, upper = 1> prior_only;
    array[N] int<lower = 1, upper = 16> pot_idx;
    array[N] int<lower = 1, upper = 3> treatment;
    array[N] int<lower = 1, upper = 3> day_idx;
    array[16] int<lower = 1, upper = 3> pot_treatment;
    
    // Outcomes ====================================
    vector[N] gammaRC;
    vector[N] phi_Po;
    vector[N] psi_Eo;
    vector[N] delta_Ro;
}
// Precompute the paired (treatment, day) 
transformed data{
    array[N] int cell_idx;
    for (i in 1:N) {
        cell_idx[i] = (day_idx[i] - 1) * 3 + treatment[i];
    }
}
// Parameters
parameters{
    // Per Submodel Matrix ======================
    matrix[3, 3] eta_rc;
    matrix[3, 3] eta_po;
    matrix[3, 3] eta_eo;
    matrix[3, 3] eta_ro;
    // Per Submodel Kappa per Day ===============
    vector<lower = 0>[3] rc_kappa_d;
    vector<lower = 0>[3] po_kappa_d;
    vector<lower = 0>[3] eo_kappa_d;
    vector<lower = 0>[3] ro_kappa_d;
    // Per Submodel Normalizing Variables =======
    vector[16] rc_zp_pot;
    vector[16] po_zp_pot;
    vector[16] eo_zp_pot;
    vector[16] ro_zp_pot;
    // Per Submodel Tau per Treatment ============
    vector<lower = 0>[3] rc_tau_t;
    vector<lower = 0>[3] po_tau_t;
    vector<lower = 0>[3] eo_tau_t;
    vector<lower = 0>[3] ro_tau_t;
}
// Parameter Transformation for Non-Centering
transformed parameters {
    // Pot-level non-centered offsets
    vector[16] rc_up = rc_tau_t[pot_treatment] .* rc_zp_pot;
    vector[16] po_up = po_tau_t[pot_treatment] .* po_zp_pot;
    vector[16] eo_up = eo_tau_t[pot_treatment] .* eo_zp_pot;
    vector[16] ro_up = ro_tau_t[pot_treatment] .* ro_zp_pot;

    // Cell-mean lookup via single-array gather on the flattened matrix
    vector[N] rc_mu = inv_logit(to_vector(eta_rc)[cell_idx] + rc_up[pot_idx]);
    vector[N] po_mu = inv_logit(to_vector(eta_po)[cell_idx] + po_up[pot_idx]);
    vector[N] eo_mu = inv_logit(to_vector(eta_eo)[cell_idx] + eo_up[pot_idx]);
    vector[N] ro_mu = inv_logit(to_vector(eta_ro)[cell_idx] + ro_up[pot_idx]);

    // Per-observation kappa gather 
    vector[N] kd_rc = rc_kappa_d[day_idx];
    vector[N] kd_po = po_kappa_d[day_idx];
    vector[N] kd_eo = eo_kappa_d[day_idx];
    vector[N] kd_ro = ro_kappa_d[day_idx];

    // Beta shape parameters, computed once and reused by both the
    vector[N] a_rc = rc_mu .* kd_rc;        vector[N] b_rc = (1 - rc_mu) .* kd_rc;
    vector[N] a_po = po_mu .* kd_po;        vector[N] b_po = (1 - po_mu) .* kd_po;
    vector[N] a_eo = eo_mu .* kd_eo;        vector[N] b_eo = (1 - eo_mu) .* kd_eo;
    vector[N] a_ro = ro_mu .* kd_ro;        vector[N] b_ro = (1 - ro_mu) .* kd_ro;
}

// Models
model{
    // GammaRC Priors ======================
    to_vector(eta_rc) ~ normal(-1, 0.3);
    rc_kappa_d        ~ exponential(0.001);
    rc_tau_t          ~ exponential(2);
    rc_zp_pot         ~ normal(0, 1);

    // Phi Po Priors =======================
    to_vector(eta_po) ~ normal(0.8, 0.3);
    po_kappa_d        ~ exponential(0.001);
    po_tau_t          ~ exponential(3);
    po_zp_pot         ~ normal(0, 1);

    // Psi Eo Priors ========================
    to_vector(eta_eo) ~ normal(0.7, 0.3);
    eo_kappa_d        ~ exponential(0.01);
    eo_tau_t          ~ exponential(3);
    eo_zp_pot         ~ normal(0, 1);

    // Delta Ro Priors =======================
    to_vector(eta_ro) ~ normal(0.1, 0.7);
    ro_kappa_d        ~ exponential(0.01);
    ro_tau_t          ~ exponential(2);
    ro_zp_pot         ~ normal(0, 1);

    // Submodels 
    if (prior_only == 0) {
        gammaRC  ~ beta(a_rc, b_rc);
        phi_Po   ~ beta(a_po, b_po);
        psi_Eo   ~ beta(a_eo, b_eo);
        delta_Ro ~ beta(a_ro, b_ro);
    }
}

// Additional Calculations
generated quantities {
    // 1 POPULATION CELL MEANS, probability scale ============
    matrix[3,3] mu_rc = inv_logit(eta_rc);
    matrix[3,3] mu_po = inv_logit(eta_po);
    matrix[3,3] mu_eo = inv_logit(eta_eo);
    matrix[3,3] mu_ro = inv_logit(eta_ro);

    // 2 DERIVED ESTIMANDS: RC/ABS, PIabs, PItotal per cell ==================
    matrix[3,3] RC_ABS  = mu_rc ./ (1 - mu_rc);
    matrix[3,3] po_odds = mu_po ./ (1 - mu_po);
    matrix[3,3] eo_odds = mu_eo ./ (1 - mu_eo);
    matrix[3,3] ro_odds = mu_ro ./ (1 - mu_ro);
    matrix[3,3] PIabs   = RC_ABS .* po_odds .* eo_odds;
    matrix[3,3] PItotal = PIabs .* ro_odds;

    // 3 CONTRASTS per day (EMF=row 3, Sham=row 2, Control=row 1) ===========
    vector[3] d_rc_emf_ctrl  = to_vector(eta_rc[3] - eta_rc[1]);
    vector[3] d_rc_emf_sham  = to_vector(eta_rc[3] - eta_rc[2]);
    vector[3] d_rc_sham_ctrl = to_vector(eta_rc[2] - eta_rc[1]);

    vector[3] d_po_emf_ctrl  = to_vector(eta_po[3] - eta_po[1]);
    vector[3] d_po_emf_sham  = to_vector(eta_po[3] - eta_po[2]);
    vector[3] d_po_sham_ctrl = to_vector(eta_po[2] - eta_po[1]);

    vector[3] d_eo_emf_ctrl  = to_vector(eta_eo[3] - eta_eo[1]);
    vector[3] d_eo_emf_sham  = to_vector(eta_eo[3] - eta_eo[2]);
    vector[3] d_eo_sham_ctrl = to_vector(eta_eo[2] - eta_eo[1]);

    vector[3] d_ro_emf_ctrl  = to_vector(eta_ro[3] - eta_ro[1]);
    vector[3] d_ro_emf_sham  = to_vector(eta_ro[3] - eta_ro[2]);
    vector[3] d_ro_sham_ctrl = to_vector(eta_ro[2] - eta_ro[1]);

    // Joint Submodel Posterior Derived Quantities Contrasts
    vector[3] logPIabs_emf_ctrl  = to_vector(log(PIabs[3])   - log(PIabs[1]));
    vector[3] logPIabs_emf_sham  = to_vector(log(PIabs[3])   - log(PIabs[2]));
    vector[3] logPIabs_sham_ctrl = to_vector(log(PIabs[2])   - log(PIabs[1]));

    vector[3] logPItot_emf_ctrl  = to_vector(log(PItotal[3]) - log(PItotal[1]));
    vector[3] logPItot_emf_sham  = to_vector(log(PItotal[3]) - log(PItotal[2]));
    vector[3] logPItot_sham_ctrl = to_vector(log(PItotal[2]) - log(PItotal[1]));

    vector[3] dPIabs_emf_ctrl    = to_vector(PIabs[3]   - PIabs[1]);
    vector[3] dPItot_emf_ctrl    = to_vector(PItotal[3] - PItotal[1]);

    // 4+5 PROBABILITY-OF-DIRECTION and ROPE indicators =====================
    real rope_logit_half = 0.10;
    real rope_logratio_half = 0.05;
    array[3] int pd_rc_emf_ctrl;      array[3] int pd_po_emf_ctrl;
    array[3] int pd_eo_emf_ctrl;      array[3] int pd_ro_emf_ctrl;
    array[3] int pd_PIabs_emf_ctrl;   array[3] int pd_PItot_emf_ctrl;
    array[3] int rope_rc_emf_ctrl;    array[3] int rope_po_emf_ctrl;
    array[3] int rope_eo_emf_ctrl;    array[3] int rope_ro_emf_ctrl;
    array[3] int rope_PIabs_emf_ctrl; array[3] int rope_PItot_emf_ctrl;
    for (d in 1:3) {
        pd_rc_emf_ctrl[d]    = is_positive(d_rc_emf_ctrl[d]);
        pd_po_emf_ctrl[d]    = is_positive(d_po_emf_ctrl[d]);
        pd_eo_emf_ctrl[d]    = is_positive(d_eo_emf_ctrl[d]);
        pd_ro_emf_ctrl[d]    = is_positive(d_ro_emf_ctrl[d]);
        pd_PIabs_emf_ctrl[d] = is_positive(logPIabs_emf_ctrl[d]);
        pd_PItot_emf_ctrl[d] = is_positive(logPItot_emf_ctrl[d]);

        rope_rc_emf_ctrl[d]    = in_rope(d_rc_emf_ctrl[d], -rope_logit_half, rope_logit_half);
        rope_po_emf_ctrl[d]    = in_rope(d_po_emf_ctrl[d], -rope_logit_half, rope_logit_half);
        rope_eo_emf_ctrl[d]    = in_rope(d_eo_emf_ctrl[d], -rope_logit_half, rope_logit_half);
        rope_ro_emf_ctrl[d]    = in_rope(d_ro_emf_ctrl[d], -rope_logit_half, rope_logit_half);
        rope_PIabs_emf_ctrl[d] = in_rope(logPIabs_emf_ctrl[d], -rope_logratio_half, rope_logratio_half);
        rope_PItot_emf_ctrl[d] = in_rope(logPItot_emf_ctrl[d], -rope_logratio_half, rope_logratio_half);
    }

    // 6 POT-LEVEL heterogeneity summaries ===================================
    real sd_rc_up = sd(rc_up); real sd_po_up = sd(po_up);
    real sd_eo_up = sd(eo_up); real sd_ro_up = sd(ro_up);

    // 7 BAYES R2 per submodel ===============================================
    vector[N] eta_rc_pop = to_vector(eta_rc)[cell_idx];
    vector[N] eta_po_pop = to_vector(eta_po)[cell_idx];
    vector[N] eta_eo_pop = to_vector(eta_eo)[cell_idx];
    vector[N] eta_ro_pop = to_vector(eta_ro)[cell_idx];
    
    // full linear predictor = population + pot offset, one vector op each.
    vector[N] eta_rc_full = eta_rc_pop + rc_up[pot_idx];
    vector[N] eta_po_full = eta_po_pop + po_up[pot_idx];
    vector[N] eta_eo_full = eta_eo_pop + eo_up[pot_idx];
    vector[N] eta_ro_full = eta_ro_pop + ro_up[pot_idx];

    // delta-method sampling variance: Var(logit) ~ 1 / ((kappa+1) mu (1-mu))
    vector[N] vrc = 1.0 ./ ((kd_rc + 1) .* rc_mu .* (1 - rc_mu));
    vector[N] vpo = 1.0 ./ ((kd_po + 1) .* po_mu .* (1 - po_mu));
    vector[N] veo = 1.0 ./ ((kd_eo + 1) .* eo_mu .* (1 - eo_mu));
    vector[N] vro = 1.0 ./ ((kd_ro + 1) .* ro_mu .* (1 - ro_mu));

    // R Squared 
    real R2_rc_full = bayes_R2_general(eta_rc_full, mean(vrc));
    real R2_rc_pop  = bayes_R2_general(eta_rc_pop,  mean(vrc) + variance(eta_rc_full - eta_rc_pop));
    real R2_po_full = bayes_R2_general(eta_po_full, mean(vpo));
    real R2_po_pop  = bayes_R2_general(eta_po_pop,  mean(vpo) + variance(eta_po_full - eta_po_pop));
    real R2_eo_full = bayes_R2_general(eta_eo_full, mean(veo));
    real R2_eo_pop  = bayes_R2_general(eta_eo_pop,  mean(veo) + variance(eta_eo_full - eta_eo_pop));
    real R2_ro_full = bayes_R2_general(eta_ro_full, mean(vro));
    real R2_ro_pop  = bayes_R2_general(eta_ro_pop,  mean(vro) + variance(eta_ro_full - eta_ro_pop));

    // 8 LOG-LIK ==============================================================
    vector[N] log_lik_rc = (a_rc - 1) .* log(gammaRC)  + (b_rc - 1) .* log1m(gammaRC)
                            - lgamma(a_rc) - lgamma(b_rc) + lgamma(a_rc + b_rc);
    vector[N] log_lik_po = (a_po - 1) .* log(phi_Po)   + (b_po - 1) .* log1m(phi_Po)
                            - lgamma(a_po) - lgamma(b_po) + lgamma(a_po + b_po);
    vector[N] log_lik_eo = (a_eo - 1) .* log(psi_Eo)   + (b_eo - 1) .* log1m(psi_Eo)
                            - lgamma(a_eo) - lgamma(b_eo) + lgamma(a_eo + b_eo);
    vector[N] log_lik_ro = (a_ro - 1) .* log(delta_Ro) + (b_ro - 1) .* log1m(delta_Ro)
                            - lgamma(a_ro) - lgamma(b_ro) + lgamma(a_ro + b_ro);
    vector[N] log_lik = log_lik_rc + log_lik_po + log_lik_eo + log_lik_ro;

    // 9 POSTERIOR PREDICTIVE REPS ===========================================
    // Small epsilon to keep mu away from exact 0/1 boundaries during prior/posterior predictive draws
    real eps = 1e-9;
    array[N] real rc_rep;       array[N] real po_rep; 
    array[N] real eo_rep;       array[N] real ro_rep;
    array[N] real PIabs_rep;    array[N] real PItot_rep;
    array[N] real logPIabs_rep; array[N] real logPItot_rep;

    for (i in 1:N) {
        real rr =  beta_rng(a_rc[i], b_rc[i]);
        real rp =  beta_rng(a_po[i], b_po[i]);
        real re =  beta_rng(a_eo[i], b_eo[i]);
        real rro = beta_rng(a_ro[i], b_ro[i]);

        // Clip each replicated draw away from the boundary before any downstream transform
        real rr_c  = fmin(fmax(rr,  eps), 1 - eps);
        real rp_c  = fmin(fmax(rp,  eps), 1 - eps);
        real re_c  = fmin(fmax(re,  eps), 1 - eps);
        real rro_c = fmin(fmax(rro, eps), 1 - eps);
        rc_rep[i] = rr; po_rep[i] = rp; eo_rep[i] = re; ro_rep[i] = rro; 

        real pia = (rr_c / (1 - rr_c)) * (rp_c / (1 - rp_c)) * (re_c / (1 - re_c));
        PIabs_rep[i] = pia;
        PItot_rep[i] = pia * (rro_c / (1 - rro_c));
        logPIabs_rep[i] = log(pia);
        logPItot_rep[i] = log(PItot_rep[i]);
    }
}
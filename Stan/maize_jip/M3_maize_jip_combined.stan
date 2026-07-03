// Function Block 
functions{
    #include ../lib/effect_probabilities.stanfunctions
    #include ../lib/diagnostics.stanfunctions
}
// Input Data 
data{
    int<lower = 0 > N;
    int<lower = 0, upper = 1 > prior_only;
    array[N] int<lower = 1, upper = 16> pot_idx;
    array[N] int<lower = 1, upper = 3> treatment;
    array[N] int<lower = 1, upper = 3> day_idx;
    array[16]int<lower = 1, upper = 3> pot_treatment;
    // Outcomes ==================================== 
    vector[N] gammaRC;
    vector[N] phi_Po;
    vector[N] psi_Eo;
    vector[N] delta_Ro;
}
// Paramteters 
parameters{
    // Per Submodel Matrix ======================
    matrix[3, 3] eta_rc;
    matrix[3, 3] eta_po;
    matrix[3, 3] eta_eo;
    matrix[3, 3] eta_ro;
    // Per Submodel Kappa per Day ===============
    vector[3] rc_kappa_d;
    vector[3] po_kappa_d;
    vector[3] eo_kappa_d;
    vector[3] ro_kappa_d;
    // Per Submodel Normalizing Variables =======
    vector[16] rc_zp_pot;
    vector[16] po_zp_pot;
    vector[16] eo_zp_pot;
    vector[16] ro_zp_pot;
    // Per Submodel Tau per Treatment ============
    vector[3] rc_tau_t;
    vector[3] po_tau_t;
    vector[3] eo_tau_t;
    vector[3] ro_tau_t;
    
}
// Parameter Transfomration for Non-Centering 
transformed parameters {
   // Caclulate U per Submodel ====================
   vector[16] rc_up;
   vector[16] po_up;
   vector[16] eo_up;
   vector[16] ro_up;
   for(p in 1:16){
    rc_up[p] = rc_tau_t[pot_treatment[p]] * rc_zp_pot[p];
    po_up[p] = po_tau_t[pot_treatment[p]] * po_zp_pot[p];
    eo_up[p] = eo_tau_t[pot_treatment[p]] * eo_zp_pot[p];
    ro_up[p] = ro_tau_t[pot_treatment[p]] * ro_zp_pot[p];
   }
   // Calculate Mu for Every submodel 
   vector[N] rc_mu;
   vector[N] po_mu;
   vector[N] eo_mu;
   vector[N] ro_mu;
   for(i in 1:N){
    rc_mu[i] = inv_logit(eta_rc[treatment[i], day_idx[i]] + rc_up[pot_idx[i]]);
    po_mu[i] = inv_logit(eta_po[treatment[i], day_idx[i]] + po_up[pot_idx[i]]);
    eo_mu[i] = inv_logit(eta_eo[treatment[i], day_idx[i]] + eo_up[pot_idx[i]]);
    ro_mu[i] = inv_logit(eta_ro[treatment[i], day_idx[i]] + ro_up[pot_idx[i]]);
   }
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
    if(prior_only == 0){
        for(i in 1:N){
            gammaRC[i]  ~ beta(rc_mu[i] * rc_kappa_d[day_idx[i]], (1 - rc_mu[i]) * rc_kappa_d[i]);
            phi_Po[i]   ~ beta(po_mu[i] * po_kappa_d[day_idx[i]], (1 - po_mu[i]) * po_kappa_d[i]);
            psi_Eo[i]   ~ beta(eo_mu[i] * eo_kappa_d[day_idx[i]], (1 - eo_mu[i]) * eo_kappa_d[i]);
            delta_Ro[i] ~ beta(ro_mu[i] * ro_kappa_d[day_idx[i]], (1 - ro_mu[i]) * ro_kappa_d[i]);
        }
    }
}
// Additional Calculations 
generated quantities {
    // 1 POPULATION CELL MEANS, probability scale  ============
    matrix[3,3] mu_rc; matrix[3,3] mu_po; matrix[3,3] mu_eo; matrix[3,3] mu_ro;
    for (t in 1:3) for (d in 1:3) {
        mu_rc[t,d] = inv_logit(eta_rc[t,d]);
        mu_po[t,d] = inv_logit(eta_po[t,d]);
        mu_eo[t,d] = inv_logit(eta_eo[t,d]);
        mu_ro[t,d] = inv_logit(eta_ro[t,d]);
    }
    // 2  DERIVED ESTIMANDS: RC/ABS PIabs PItotal per cell ====================
    //    RC/ABS  = mu_rc/(1-mu_rc)
    //    PIabs   = RC/ABS * [po/(1-po)] * [eo/(1-eo)]
    //    PItotal = PIabs  * [ro/(1-ro)]
    matrix[3,3] RC_ABS; matrix[3,3] PIabs; matrix[3,3] PItotal;
    for (t in 1:3) for (d in 1:3) {
        RC_ABS[t,d]  = mu_rc[t,d] / (1 - mu_rc[t,d]);
        real po_odds = mu_po[t,d] / (1 - mu_po[t,d]);
        real eo_odds = mu_eo[t,d] / (1 - mu_eo[t,d]);
        real ro_odds = mu_ro[t,d] / (1 - mu_ro[t,d]);
        PIabs[t,d]   = RC_ABS[t,d] * po_odds * eo_odds;
        PItotal[t,d] = PIabs[t,d] * ro_odds;
    }
    // 3 CONTRASTS per day (EMF=3, Sham=2, Control=1) =======================
    vector[3] d_rc_emf_ctrl; vector[3] d_rc_emf_sham; vector[3] d_rc_sham_ctrl;
    vector[3] d_po_emf_ctrl; vector[3] d_po_emf_sham; vector[3] d_po_sham_ctrl;
    vector[3] d_eo_emf_ctrl; vector[3] d_eo_emf_sham; vector[3] d_eo_sham_ctrl;
    vector[3] d_ro_emf_ctrl; vector[3] d_ro_emf_sham; vector[3] d_ro_sham_ctrl;
    for (d in 1:3) {
        d_rc_emf_ctrl[d]=eta_rc[3,d]-eta_rc[1,d]; d_rc_emf_sham[d]=eta_rc[3,d]-eta_rc[2,d]; d_rc_sham_ctrl[d]=eta_rc[2,d]-eta_rc[1,d];
        d_po_emf_ctrl[d]=eta_po[3,d]-eta_po[1,d]; d_po_emf_sham[d]=eta_po[3,d]-eta_po[2,d]; d_po_sham_ctrl[d]=eta_po[2,d]-eta_po[1,d];
        d_eo_emf_ctrl[d]=eta_eo[3,d]-eta_eo[1,d]; d_eo_emf_sham[d]=eta_eo[3,d]-eta_eo[2,d]; d_eo_sham_ctrl[d]=eta_eo[2,d]-eta_eo[1,d];
        d_ro_emf_ctrl[d]=eta_ro[3,d]-eta_ro[1,d]; d_ro_emf_sham[d]=eta_ro[3,d]-eta_ro[2,d]; d_ro_sham_ctrl[d]=eta_ro[2,d]-eta_ro[1,d];
    }
    vector[3] logPIabs_emf_ctrl; vector[3] logPIabs_emf_sham; vector[3] logPIabs_sham_ctrl;
    vector[3] logPItot_emf_ctrl; vector[3] logPItot_emf_sham; vector[3] logPItot_sham_ctrl;
    vector[3] dPIabs_emf_ctrl;   vector[3] dPItot_emf_ctrl;
    for (d in 1:3) {
        logPIabs_emf_ctrl[d]  = log(PIabs[3,d])   - log(PIabs[1,d]);
        logPIabs_emf_sham[d]  = log(PIabs[3,d])   - log(PIabs[2,d]);
        logPIabs_sham_ctrl[d] = log(PIabs[2,d])   - log(PIabs[1,d]);
        logPItot_emf_ctrl[d]  = log(PItotal[3,d]) - log(PItotal[1,d]);
        logPItot_emf_sham[d]  = log(PItotal[3,d]) - log(PItotal[2,d]);
        logPItot_sham_ctrl[d] = log(PItotal[2,d]) - log(PItotal[1,d]);
        dPIabs_emf_ctrl[d]    = PIabs[3,d]   - PIabs[1,d];
        dPItot_emf_ctrl[d]    = PItotal[3,d] - PItotal[1,d];
    }
    // 4 PROBABILITY-OF-DIRECTION indicators =============================
    array[3] int pd_rc_emf_ctrl; array[3] int pd_po_emf_ctrl;
    array[3] int pd_eo_emf_ctrl; array[3] int pd_ro_emf_ctrl;
    array[3] int pd_PIabs_emf_ctrl; array[3] int pd_PItot_emf_ctrl;
    for (d in 1:3) {
        pd_rc_emf_ctrl[d]    = is_positive(d_rc_emf_ctrl[d]);
        pd_po_emf_ctrl[d]    = is_positive(d_po_emf_ctrl[d]);
        pd_eo_emf_ctrl[d]    = is_positive(d_eo_emf_ctrl[d]);
        pd_ro_emf_ctrl[d]    = is_positive(d_ro_emf_ctrl[d]);
        pd_PIabs_emf_ctrl[d] = is_positive(logPIabs_emf_ctrl[d]);
        pd_PItot_emf_ctrl[d] = is_positive(logPItot_emf_ctrl[d]);
    }
    // 5 ROPE indicators =================================================
    real rope_logit_half = 0.10;    // marginal primaries (logit scale)
    real rope_logratio_half = 0.05; // derived indices (log-ratio scale)
    array[3] int rope_rc_emf_ctrl; array[3] int rope_po_emf_ctrl;
    array[3] int rope_eo_emf_ctrl; array[3] int rope_ro_emf_ctrl;
    array[3] int rope_PIabs_emf_ctrl; array[3] int rope_PItot_emf_ctrl;
    for (d in 1:3) {
        rope_rc_emf_ctrl[d]    = in_rope(d_rc_emf_ctrl[d], -rope_logit_half, rope_logit_half);
        rope_po_emf_ctrl[d]    = in_rope(d_po_emf_ctrl[d], -rope_logit_half, rope_logit_half);
        rope_eo_emf_ctrl[d]    = in_rope(d_eo_emf_ctrl[d], -rope_logit_half, rope_logit_half);
        rope_ro_emf_ctrl[d]    = in_rope(d_ro_emf_ctrl[d], -rope_logit_half, rope_logit_half);
        rope_PIabs_emf_ctrl[d] = in_rope(logPIabs_emf_ctrl[d], -rope_logratio_half, rope_logratio_half);
        rope_PItot_emf_ctrl[d] = in_rope(logPItot_emf_ctrl[d], -rope_logratio_half, rope_logratio_half);
    }
    // 6 POT-LEVEL heterogeneity summaries ================================ 
    real sd_rc_up = sd(rc_up); real sd_po_up = sd(po_up);
    real sd_eo_up = sd(eo_up); real sd_ro_up = sd(ro_up);

    // 7 BAYES R2 per submodel ==============================================
    vector[N] eta_rc_full; vector[N] eta_rc_pop;
    vector[N] eta_po_full; vector[N] eta_po_pop;
    vector[N] eta_eo_full; vector[N] eta_eo_pop;
    vector[N] eta_ro_full; vector[N] eta_ro_pop;
    vector[N] vrc; vector[N] vpo; vector[N] veo; vector[N] vro;
    for (i in 1:N) {
        int d = day_idx[i]; int t = treatment[i];
        eta_rc_full[i]=eta_rc[t,d]+rc_up[pot_idx[i]]; eta_rc_pop[i]=eta_rc[t,d];
        eta_po_full[i]=eta_po[t,d]+po_up[pot_idx[i]]; eta_po_pop[i]=eta_po[t,d];
        eta_eo_full[i]=eta_eo[t,d]+eo_up[pot_idx[i]]; eta_eo_pop[i]=eta_eo[t,d];
        eta_ro_full[i]=eta_ro[t,d]+ro_up[pot_idx[i]]; eta_ro_pop[i]=eta_ro[t,d];
        // delta-method: Var(logit) ~ Var_beta(mu) / (mu(1-mu))^2 = 1/((k+1) mu (1-mu))
        real mrc=inv_logit(eta_rc_full[i]); vrc[i]=1.0/((rc_kappa_d[d]+1)*mrc*(1-mrc));
        real mpo=inv_logit(eta_po_full[i]); vpo[i]=1.0/((po_kappa_d[d]+1)*mpo*(1-mpo));
        real meo=inv_logit(eta_eo_full[i]); veo[i]=1.0/((eo_kappa_d[d]+1)*meo*(1-meo));
        real mro=inv_logit(eta_ro_full[i]); vro[i]=1.0/((ro_kappa_d[d]+1)*mro*(1-mro));
    }
    real R2_rc_full = bayes_R2_general(eta_rc_full, mean(vrc));
    real R2_rc_pop  = bayes_R2_general(eta_rc_pop,  mean(vrc)+variance(eta_rc_full-eta_rc_pop));
    real R2_po_full = bayes_R2_general(eta_po_full, mean(vpo));
    real R2_po_pop  = bayes_R2_general(eta_po_pop,  mean(vpo)+variance(eta_po_full-eta_po_pop));
    real R2_eo_full = bayes_R2_general(eta_eo_full, mean(veo));
    real R2_eo_pop  = bayes_R2_general(eta_eo_pop,  mean(veo)+variance(eta_eo_full-eta_eo_pop));
    real R2_ro_full = bayes_R2_general(eta_ro_full, mean(vro));
    real R2_ro_pop  = bayes_R2_general(eta_ro_pop,  mean(vro)+variance(eta_ro_full-eta_ro_pop));

    // 8. LOG-LIK and REPS ===================================================================
    vector[N] log_lik;
    vector[N] log_lik_rc; 
    vector[N] log_lik_po; 
    vector[N] log_lik_eo; 
    vector[N] log_lik_ro;
    array[N] real rc_rep;
    array[N] real po_rep; 
    array[N] real eo_rep; 
    array[N] real ro_rep;
    array[N] real PIabs_rep; 
    array[N] real PItot_rep;

    for (i in 1:N) {
        int d = day_idx[i];
        real arc = rc_mu[i] * rc_kappa_d[d]; 
        real brc = (1 - rc_mu[i]) * rc_kappa_d[d];
        real apo = po_mu[i] * po_kappa_d[d];
        real bpo = (1 - po_mu[i]) * po_kappa_d[d];
        real aeo = eo_mu[i] * eo_kappa_d[d]; 
        real beo =(1 - eo_mu[i]) * eo_kappa_d[d];
        real aro = ro_mu[i] * ro_kappa_d[d];
        real bro = (1 - ro_mu[i]) * ro_kappa_d[d];
        log_lik_rc[i] = beta_lpdf(gammaRC[i] |arc,brc);
        log_lik_po[i] = beta_lpdf(phi_Po[i]  |apo,bpo);
        log_lik_eo[i] = beta_lpdf(psi_Eo[i]  |aeo,beo);
        log_lik_ro[i] = beta_lpdf(delta_Ro[i]|aro,bro);
        log_lik[i]    = log_lik_rc[i] + log_lik_po[i] + log_lik_eo[i] + log_lik_ro[i];
        real rr = beta_rng(arc,brc);
        real rp = beta_rng(apo,bpo);
        real re = beta_rng(aeo,beo);
        real ro= beta_rng(aro,bro);
        rc_rep[i] = rr; 
        po_rep[i] = rp; 
        eo_rep[i] = re;
        ro_rep[i] = ro;
        real pia = (rr/(1-rr))*(rp/(1-rp))*(re/(1-re));
        PIabs_rep[i] = pia;
        PItot_rep[i] = pia * (ro / (1 - ro));
    }
}
/* 
M3 -- adds beta_sod_dw allowed to vary by day (SOD -> H2O2), on top of M2
 SOD -> H2O2 is the enzymatic step: SOD converts superoxide into H2O2 as
 part of the plant's active defense. This is a REGULATED response, not a
 passive chemical process, so its strength can plausibly change as the
 defense program itself changes over time:

 - Early: SOD activity may be ramping up in response to a new stressor
    -> a strong, fast-responding SOD->H2O2 coupling.
 - Later: SOD activity may plateau, be down-regulated (defense
    already established / stressor no longer novel), or become
    decoupled from H2O2 if other pathways (CAT, TEAC) take over primary
    detoxification duty -> potentially a weaker or different coupling.
    CAT->H2O2 and TEAC->H2O2 are deliberately left FIXED (not day-varying)
    in both M2 and M3: their real-data correlations with H2O2 are already
    weak, and adding day-variation to an already-weak, already-uncertain
    edge would ask more of the data (n=48, 16 pots) than it can support.
*/
// Function Block
functions{
    #include ../lib/effect_probabilities.stanfunctions
    #include ../lib/diagnostics.stanfunctions
    #include ../lib/utilities.stanfunctions
}
// DATA INPUT BLOCK : Settings + Outcomes and Covariates 
data{
    // SETTINGS ========================================
    /*
    This is CROSSED UNBALANCED DESIGN so the indexes are:
    1. cell_id is the identifier for the combination matrix of unique Treatment x Day
    2. pot_treatment_id is the identifier for the combination of unique Pot x Treatment
    */
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;
    array[N]  int <lower = 1, upper = 3>   treatment_id;
    array[N]  int <lower = 1, upper = 3>   day_id;
    array[N]  int <lower = 1, upper = 9>   cell_id;
    array[N]  int <lower = 1, upper = 16>  pot_id;     
    array[16] int <lower = 1, upper = 16>  pot_treatment_id;  

    // OUTCOMES PRIMARY BLOCK ===========================
    /*
    All of the primary block outcomes are estimated individual models and they dont relate to each other
    This block contains SOD, CAT, TEAC, Sugars.
    Every one of them we assume Treatment    -> Outcome (Direct Treatment effect)
    None of them can be negative values so   -> lower = 0
    And all of are measures of Dry Weight so -> DW
    */ 
    vector <lower = 0>[N] sod_dw;
    vector <lower = 0>[N] cat_dw;  
    vector <lower = 0>[N] sugar_dw;
    vector <lower = 0>[N] trolox_dw;

    // OUTCOMES SECONDARY BLOCK =========================
    /*
    In this block we have only H202. We assume the following :
    Treatment -> SOD  -> H202
    Treatment -> CAT  -> H202
    Treatment -> TEAC -> H202
    So the treatment effect is indirect and is passing from the primaries above.
    */ 
    vector <lower = 0>[N] h2o2_dw;

    // OUTCOMES TERTIARY BLOCK ==========================
    /*
    In this block we only have MDA. Again the Assumption is :
    H202 -> MDA
    So the treatment effect passes from each of the primaries into the H202 then into MDA
    */
    vector <lower = 0>[N] mda_dw;

    // OUTCOME BIOMASS BLOCK ============================
    /*
    In this block we have the final estimand of interest which is : Grams of water present, per gram of dry (solid) tissue
    Assumption is :
    Sugar -> Water Content ... where the Sugar is a primary such as : Treatment -> Sugar
    MDA   -> Water Content  
    */
    vector <lower = 0>[N] water_content_dw;
}
// TRANSFORMED DATA 
transformed data {
   // Get the mean and sd from the data computed once and use them later in the scale_fix 
   array[7] real samples_log_mean;
   array[7] real samples_log_sd;
   {
    // Combine the data in one array to loop over
    array[7] vector[N] all_samples = {sod_dw, cat_dw, sugar_dw, trolox_dw, h2o2_dw, mda_dw, water_content_dw}; 
    
    for(i in 1:7){
        samples_log_mean[i] = mean(log(all_samples[i]));
        samples_log_sd[i]   = sd(log(all_samples[i]));
    }
   }
}
// MODELS PARAMETERS 
parameters{
    // RANDOM INTERCEPT PARAMETERS ==================================================================================
    // All of the parameters below are required for construction of the random partially pooled intercepts as:
    // mu = ni[t(i), d(i)] + u[p]... This model structure is common for all Submodels

    // PARAMETERS Cell Matrix =========================
    // Every Submodel gets its own Cell Matrix 3 x 3. It for estimating 9 intercepts one per each Treatment x Day Combination cell
    array[7] matrix[3, 3] eta; 

    // PARAMETERS ZETA_P Normalizing parameter =========
    // Every Submodel Gets its own Zeta parameters required for Non Centered Parametrization (NCP).Canonically zeta_x ~ N(0, 1)
    // The Zeta Normalization lives in the pot level so vector of length 16
    array[7] vector[16] zeta_p;

    // PARAMETERS TAU_T Variation ======================
    // Every Submodel gets it own Tau_t variation that allows each treatment to have its own variance
    // The Tau_t lives in the treatment level in our case 3 treatments Control SHAM EMF so Tau_t is vector of length 3
    array[7] vector <lower = 0>[3] tau_t;    

    // PARAMETERS SD_D Observation ======================
    // Every Submodel gets it own SD_D variation that allows per day variance
    // The SD_D lives in the day level in our we have 3 measurements day 14 20 28 so SD_D is a vector of length 3
    array[7] vector <lower = 0>[3] sd_obs;

    // RELATIONSHIP PARAMETERS ====================================================================================
    // This parameters encodes the causal relations between the models discussed the data input block
    
    // CAT -> H202              TEAC -> H202
    real beta_cat_dw;           real beta_trolox_dw;
    
    // Sugar -> Water Content   MDA -> Water Content  
    real beta_sugar_dw;         real beta_mda_dw;
    
    // Random Effect H202 -> MDA NCP with domain that lives in the day    
    real              beta_mean_h2o2_dw;
    real <lower = 0>  sd_h2o2_dw;
    vector[3]         zd_h2o2_dw; 

    // Random NCP SOD -> H202, NCP with domain that lives in the day      
    real                beta_mean_sod_dw;
    real <lower = 0.01> beta_sd_sod_dw;
    vector[3]           zd_sod_dw;
}
// TRANSFORM PARAMETERS
transformed parameters {
    // Recover the random partially polled intercept for all Submodels =============================================
    // Pot-level non-centered offsets
    array[7] vector[16] pot_offsets;
    array[7] vector[N]  ni;
    for(i in 1:7){
        pot_offsets[i] = nco(tau_t[i],  pot_treatment_id, zeta_p[i]);
        ni[i]          = get_ni(eta[i], cell_id, pot_id, pot_offsets[i]);
    }

    // Recover the random Effect on Effect H202 -> MDA , NCP
    vector[3] beta_h2o2_dw = beta_mean_h2o2_dw  + sd_h2o2_dw      * zd_h2o2_dw;
    
    // Recover  SOD -> H202, NCP     
    vector[3] beta_sod_dw  = beta_mean_sod_dw   + beta_sd_sod_dw  * zd_sod_dw;

    // Scale every submodel ni that will go in another model mu and for the primaries ni is mu
    array[4] vector[N] primaries_mu_z;
    {
      for(i in 1:4){
        primaries_mu_z[i] = scale_fixed(ni[i], samples_log_mean[i], samples_log_sd[i]);  
      }
    }
    // Construct the rest of the linear predictors
    vector[N] meanlog_h2o2 =
        ni[5]
        + beta_sod_dw[day_id]     .* primaries_mu_z[1]
        + beta_cat_dw              * primaries_mu_z[2]
        + beta_trolox_dw           * primaries_mu_z[4];

    vector[N] meanlog_h2o2_z = scale_fixed(meanlog_h2o2, samples_log_mean[5], samples_log_sd[5]);

    vector[N] meanlog_mda = ni[6] + beta_h2o2_dw[day_id] .* meanlog_h2o2_z;

    vector[N] meanlog_mda_z = scale_fixed(meanlog_mda, samples_log_mean[6], samples_log_sd[6]);

    vector[N] meanlog_water =
        ni[7]
        + beta_mda_dw   * meanlog_mda_z
        + beta_sugar_dw * primaries_mu_z[3];
}
// Model Priors and Likelihood 
model{

    // PRIORS ======================================================================================

    // PRIORS Cell Matrix =========================
    to_vector(eta[1]) ~ normal(5.31, 0.25); to_vector(eta[5]) ~ normal(4.89, 0.27); 
    to_vector(eta[2]) ~ normal(0.79, 0.17); to_vector(eta[6]) ~ normal(5.78, 0.13);
    to_vector(eta[3]) ~ normal(4.45, 0.11); to_vector(eta[7]) ~ normal(2.26, 0.09);
    to_vector(eta[4]) ~ normal(5.02, 0.41);

    // PRIORS ZETA_P Normalizing parameter =========
    for (i in 1:7) zeta_p[i] ~ normal(0, 1);

    // PRIORS TAU_T Variation ======================
    tau_t[1] ~ exponential(15); tau_t[4] ~ exponential(15); tau_t[6]  ~ exponential(30);
    tau_t[2] ~ exponential(10); tau_t[5] ~ exponential(20); tau_t[7]  ~ exponential(15);
    tau_t[3] ~ exponential(15);

    // PRIORS SD_D Observation ======================
    sd_obs[1] ~ exponential(10.5); sd_obs[2] ~ exponential(3.4); sd_obs[3] ~ exponential(3.6);
    sd_obs[4] ~ exponential(12);   sd_obs[5] ~ exponential(5.8);  sd_obs[6] ~ exponential(5.3);
    sd_obs[7] ~ exponential(11);

    // RELATIONSHIP PRIORS ====================================================================================
    // CAT -> H202                       TEAC -> H202
    beta_cat_dw ~ normal(0.05, 0.08);    beta_trolox_dw ~ normal(0.05, 0.08);
    
    // Sugar -> Water Content            MDA -> Water Content  
    beta_sugar_dw  ~ normal(0, 0.10);    beta_mda_dw  ~ normal(0, 0.10); 

    // Random Effect H202 -> MDA NCP with domain that lives in the day    
    beta_mean_h2o2_dw ~ normal(0.18, 0.10);
    sd_h2o2_dw        ~ exponential(10);
    zd_h2o2_dw        ~ normal(0, 1); 

    // Random Effect SOD -> H202  NCP with domain that lives in the day    
    beta_mean_sod_dw ~ normal(0.21, 0.10);
    beta_sd_sod_dw   ~ exponential(10);
    zd_sod_dw        ~ normal(0, 1);

    // Models Likelihoods
    if(prior_only == 0){
        // PRIMARY MODELS =====================================================
        sod_dw    ~ lognormal(ni[1], sd_obs[1][day_id]);
        cat_dw    ~ lognormal(ni[2], sd_obs[2][day_id]);
        sugar_dw  ~ lognormal(ni[3], sd_obs[3][day_id]);
        trolox_dw ~ lognormal(ni[4], sd_obs[4][day_id]);

        // SECONDARY MODEL ====================================================
        h2o2_dw ~ lognormal(meanlog_h2o2, sd_obs[5][day_id]);

        // TERTIARY MODEL =====================================================
        mda_dw ~ lognormal(meanlog_mda, sd_obs[6][day_id]);

        // BIOMASS MODEL ======================================================
        water_content_dw ~ lognormal(meanlog_water, sd_obs[7][day_id]);
    }
}
generated quantities {
    array[7] vector[N] meanlog;
    array[7] vector[N] samples_rep;
    array[7] vector[N] log_lik;
    array[7] vector[N] all_samples = {sod_dw, cat_dw, sugar_dw, trolox_dw, h2o2_dw, mda_dw, water_content_dw}; 
    for (n in 1:N) {
        for (k in 1:4) {
            meanlog[k][n] = ni[k][n];
        }
        meanlog[5][n] = meanlog_h2o2[n];
        meanlog[6][n] = meanlog_mda[n];
        meanlog[7][n] = meanlog_water[n];
        
        for (k in 1:7) {
            samples_rep[k][n] = lognormal_rng(meanlog[k][n], sd_obs[k][day_id[n]]);
            log_lik[k][n]     = lognormal_lpdf(all_samples[k][n] | meanlog[k][n], sd_obs[k][day_id[n]]);
        }
    }
}
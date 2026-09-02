/* 
  The base model (M1) assumes each mediation path (e.g. SOD -> H2O2,
   H2O2 -> MDA) has a SINGLE fixed strength across the whole 14/21/28-day
   window. That's a simplifying assumption, not a biological fact -- these
   are oxidative-stress pathways, and their coupling strength plausibly
   shifts over time as the plant's antioxidant response develops. M2/M3
   test whether the data actually supports that added flexibility, rather
   than assuming it a priori.
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

   // Get the 25th and 75th percentile of every covariant in the data
   int idx_25 = to_int(round(0.25 * N));
   int idx_75 = to_int(round(0.75 * N));

   array[7] real all_samples_z_q25;
   array[7] real all_samples_z_q75;
   {
    // Combine the data in one array to loop over
    array[7] vector[N] all_samples = {sod_dw, cat_dw, sugar_dw, trolox_dw, h2o2_dw, mda_dw, water_content_dw}; 

    // Store intermediate results 
    array[7] vector[N] all_samples_sorted;
    array[7] real      all_samples_q25;
    array[7] real      all_samples_q75;

    for(i in 1:7){
        // Calculate the log mean and log sd 
        samples_log_mean[i] = mean(log(all_samples[i]));
        samples_log_sd[i]   = sd(log(all_samples[i]));

        // Calculate the percentiles 25th and 75th Use in generated quantities
        all_samples_sorted[i] = sort_asc(log(all_samples[i]));
        all_samples_q25[i]    = all_samples_sorted[i][idx_25];
        all_samples_q75[i]    = all_samples_sorted[i][idx_75];

        // Z-score the q25 and q75 for all samples
        all_samples_z_q25[i]  = (all_samples_q25[i] - samples_log_mean[i]) / samples_log_sd[i];
        all_samples_z_q75[i]  = (all_samples_q75[i] - samples_log_mean[i]) / samples_log_sd[i];
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
    // SOD -> H202               CAT -> H202             TEAC -> H202
    real beta_sod_dw;           real beta_cat_dw;        real beta_trolox_dw;
    // H202 -> MDA              Sugar -> Water Content   MDA -> Water Content  
    real beta_h2o2_dw;          real beta_sugar_dw;      real beta_mda_dw;
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
        + beta_sod_dw    * primaries_mu_z[1]
        + beta_cat_dw    * primaries_mu_z[2]
        + beta_trolox_dw * primaries_mu_z[4];

    vector[N] meanlog_h2o2_z = scale_fixed(meanlog_h2o2, samples_log_mean[5], samples_log_sd[5]);

    vector[N] meanlog_mda = ni[6] + beta_h2o2_dw * meanlog_h2o2_z;

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
    // SOD -> H202                       CAT -> H202                       TEAC -> H202
    beta_sod_dw    ~ normal(0.21, 0.10); beta_cat_dw ~ normal(0.05, 0.08); beta_trolox_dw ~ normal(0.05, 0.08); 
    
    // H202 -> MDA                      Sugar -> Water Content             MDA -> Water Content  
    beta_h2o2_dw   ~ normal(0.18, 0.10); beta_sugar_dw  ~ normal(0, 0.10); beta_mda_dw    ~ normal(0, 0.10); 

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
    // POPULATION CELL MEANS, outcome scale =================================================================
    array[7] vector[N] mu_outcome_scale;  
    for(i in 1:7){
       mu_outcome_scale[i] = exp(to_vector(eta[i]));
    } 

    // CONTRASTS per day (EMF=row 3, Sham=row 2, Control=row 1)==============================================
    array[7] vector[3] d_emf_ctrl;
    array[7] vector[3] d_emf_sham;
    array[7] vector[3] d_sham_ctrl;
    for(i in 1:7){
        d_emf_ctrl[i]  = to_vector(eta[i][3] - eta[i][1]); 
        d_emf_sham[i]  = to_vector(eta[i][3] - eta[i][2]);
        d_sham_ctrl[i] = to_vector(eta[i][2] - eta[i][1]);
    }

    // DERIVED QUANTITIES ===================================================================================

    // SOD-associated H₂O₂ production pressure
    //  How strongly the SOD pathway pushes the system toward H₂O₂ under the specified SOD change.
    real delta_H_sod_dw = beta_sod_dw * (all_samples_z_q75[1] - all_samples_z_q25[1]);
    
    // CAT-mediated H₂O₂ detoxification
    //  How much H₂O₂ reduction is associated with the specified increase in CAT.
    real delta_H_cat_dw = beta_cat_dw * (all_samples_z_q75[2] - all_samples_z_q25[2]);
    
    // TEAC-mediated H₂O₂ protection
    // How strongly the non-enzymatic antioxidant capacity represented by TEAC suppresses H₂O₂ under the specified TEAC change
    real delta_H_teac_dw = beta_trolox_dw * (all_samples_z_q75[4] - all_samples_z_q25[4]);
    
    // ROS buffering ratio
    // Relative CAT buffering strength compared with SOD-associated H₂O₂ pressure.
    real ros_buffer_ration = delta_H_cat_dw / delta_H_teac_dw;

    // Relative antioxidant contribution
    // Of the combined modeled CAT + TEAC H₂O₂-reducing response, how much is attributable to TEAC
    real R_antiox = delta_H_teac_dw / (delta_H_cat_dw + delta_H_teac_dw);

    // Oxidative damage susceptibility
    // H₂O₂ increases by the specified amount, how much does the model predict MDA to increase?
    real delta_H_h202_dw = beta_h2o2_dw * (all_samples_z_q75[5] - all_samples_z_q25[5]); 

    // Damage-associated physiological effect
    // If lipid damage increases, how much does the model predict water content to change?
    real delta_biomass_mda = beta_mda_dw * (all_samples_z_q75[6] - all_samples_z_q25[6]);

    // MEDIATION ============================================================================================
    /* Downstream   Mediator  "Direct" bucket              	Question this answers
1	     H₂O₂	    SOD	       CAT + TEAC + H₂O₂'s own eta	How much of Treatment's effect on H₂O₂ flows through SOD specifically?
2	     H₂O₂	    CAT        SOD + TEAC + H₂O₂'s own eta	How much flows through CAT specifically?
3	     H₂O₂	    TEAC	   SOD + CAT + H₂O₂'s own eta	How much flows through TEAC specifically?
4	     MDA	    H₂O₂	   MDA's own eta	            How much of Treatment's effect on MDA flows through H₂O₂?
5	     Water	    MDA	       Sugars + Water's own eta	    How much of Treatment's effect on Water flows through MDA?
6	     Water	    Sugars	   MDA + Water's own eta	    How much of Treatment's effect on Water flows through Sugars?
    */ 
    // Estimands of interest 
    array[3, 3, 3] real PNDE_h2o2;  array[3, 3, 3] real PNDE_mda; array[3, 3, 3] real PNDE_water;
    array[3, 3, 3] real TNIE_h2o2;  array[3, 3, 3] real TNIE_mda; array[3, 3, 3] real TNIE_water;
    
    {
        // Declare Z-score matrixes for every cell model
        array[4] matrix[3, 3] primaries_cell_z;
        matrix[3, 3] mda_cell;    matrix[3, 3] mda_cell_z; 
        matrix[3, 3] h2_o2_cell;  matrix[3, 3] h2_o2_cell_z;
        matrix[3, 3] water_cell;  matrix[3, 3] water_cell_z;

        // Loop over every treatment and day to fill the matrixes
        for(t in 1:3){
            for(d in 1:3){
                // Zscore the primaries
                for(p in 1:4){
                    primaries_cell_z[t, d, p] = (eta[p][t, d] - samples_log_mean[p]) / samples_log_sd[p];
                } 
                // Compute the h2o2 cell matrix with the Z-scored primary cell matrixes and Z-score the h2o2 matrix
                h2_o2_cell[t, d] = eta[5][t, d] 
                                    + beta_sod_dw    * primaries_cell_z[1][t, d]
                                    + beta_cat_dw    * primaries_cell_z[2][t, d] 
                                    + beta_trolox_dw * primaries_cell_z[4][t, d];
                
                h2_o2_cell_z[t, d]  = (h2_o2_cell[t, d] - samples_log_mean[5]) / samples_log_sd[5];

                // Compute the MDA cell matrix and Z-score it 
                mda_cell[t, d]   = eta[6][t, d]    + beta_h2o2_dw         * h2_o2_cell_z[t, d];
                mda_cell_z[t, d] = (mda_cell[t, d] - samples_log_mean[6]) / samples_log_sd[6];

                // Compute water cell matrix and Z-score it 
                water_cell[t, d] = eta[7][t, d] 
                                    + beta_mda_dw   * mda_cell_z[t, d] 
                                    + beta_sugar_dw * primaries_cell_z[3][t, d];
                water_cell_z[t, d] = (water_cell[t, d] - samples_log_mean[7]) / samples_log_sd[7];
            }
        }
        // Compute PNDE and TNIE 
        for (t1 in 1:3) for (t0 in 1:3) for (d in 1:3) {
            // H2o2
            PNDE_h2o2[t1, t0, d] = eta[5][t1, d] - eta[5][t0, d];
            TNIE_h2o2[t1, t0, d] = beta_sod_dw       * (primaries_cell_z[1][t1, d] - primaries_cell_z[1][t0, d])
                                    + beta_cat_dw    * (primaries_cell_z[2][t1, d] - primaries_cell_z[2][t0, d])
                                    + beta_trolox_dw * (primaries_cell_z[4][t1, d] - primaries_cell_z[4][t0, d]);
            // MDA
            PNDE_mda[t1, t0, d] =  eta[6][t1, d] - eta[6][t0, d];
            TNIE_mda[t1, t0, d] =  beta_h2o2_dw  * (h2_o2_cell_z[t1, d] - h2_o2_cell_z[t0, d]);

            // Water 
            PNDE_water[t1, t0, d] = eta[7][t1, d]  - eta[7][t0, d];
            TNIE_water[t1, t0, d] = beta_mda_dw    * (mda_cell_z[t1, d]          - mda_cell_z[t0, d])
                                   + beta_sugar_dw * (primaries_cell_z[3][t1, d] - primaries_cell_z[3][t0, d]);
        }
    }

    // REPLICATIONS  and COMPARISON LOO =====================================================================
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
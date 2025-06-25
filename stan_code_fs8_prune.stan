// Version History. Update in print statement in transformed data!
// Version 1.5. First official release with new selection-effect and population model!
// Version 1.6 (Sep-17-2020). Added full k-correction uncertainty propagation to selection-effect model.
// Version 1.7 (Jan-26-2024). Fixed parameter limits on selection-effect model outl_mBx1c_uncertainties.
// Version 1.71 (Sep-13-2024). Added lower limit to mobs_var_by_SN_except_c_R Thanks Aaron Do!

//// MOVE CALCULATIONS TO FUNCTIONS AND INCLUDE CALCULATION OF Hr, which is only implemented for cosmo_model==1

functions{
   array[] vector calc_model_mu_Hr(int n_sne, int cosmo_model, int nzadd, real Om, array[] real redshifts, array[] real redshifts_sort_fill, array[] int unsort_inds, array[] real zhelio, array[] int photoz_inds){
        array[2] vector[n_sne] out;
        vector [n_sne]  model_mu;
        vector [n_sne]  Hr;
        array[2*(n_sne + nzadd) - 1] real Hinv_sort_fill;
        array[n_sne + nzadd] real r_com_sort;
        real dz_term;
        real dz_Hinv_term;       

        // temporary variables that should be arguments when running something other than cosmo_model=1
        real wDE;
        real waDE;
        real q0;
        real j0;
        int n_photoz=0;
        int n_zbins=0;
        vector [n_photoz] dz;
        matrix [n_sne, n_zbins] dmu_dbin;
        matrix [n_sne, n_zbins] dmudz_dbin;
        vector [n_zbins] mu_zbins;
        vector [n_sne] mu_const;
        vector [n_zbins] r_comove_bins;
        vector [n_zbins] zbins;

        // -------------Begin numerical integration-----------------


        if ((cosmo_model == 1) || (cosmo_model == 3) || (cosmo_model == 5)) {
            for (i in 1: 2*(n_sne + nzadd) - 1) {    // Inverse Hubble parameter
                if (cosmo_model == 1) {
                    Hinv_sort_fill[i] = 1./sqrt( Om*pow(1. + redshifts_sort_fill[i], 3) + (1. - Om) );
                }
                if (cosmo_model == 3) {
                    Hinv_sort_fill[i] = 1./sqrt( Om*pow(1. + redshifts_sort_fill[i], 3) + (1. - Om)*pow(1. + redshifts_sort_fill[i], 3.*(1 + wDE)) );
                }
                if (cosmo_model == 5) {
                    Hinv_sort_fill[i] = 1./sqrt( Om*pow(1. + redshifts_sort_fill[i], 3)
                        + (1. - Om)*pow(1. + redshifts_sort_fill[i], 3.*(1 + wDE + waDE))*exp(-3.*waDE*redshifts_sort_fill[i]/(1. + redshifts_sort_fill[i])) );
                }
            }

            // Integrate comoving r using Simpson's rule

            r_com_sort[1] = 0.; // Redshift = 0 should be first element!
            for (i in 2:(n_sne + nzadd)) {
                r_com_sort[i] = r_com_sort[i - 1] + (Hinv_sort_fill[2*i - 3] + 4.*Hinv_sort_fill[2*i - 2] + Hinv_sort_fill[2*i - 1])*(redshifts_sort_fill[2*i - 1] - redshifts_sort_fill[2*i - 3])/6.;
            }


            for (i in 1:n_sne) {
                if (photoz_inds[i] == 0) {
                    dz_term = 0.;
                    dz_Hinv_term = 0.;
                } 
                else {
                    dz_term = dz[photoz_inds[i]];
                    if (cosmo_model == 1) {
                        dz_Hinv_term = dz[photoz_inds[i]]/sqrt( Om*pow(1. + redshifts[i], 3) + (1. - Om) );
                    }
                    if (cosmo_model == 3) {
                        dz_Hinv_term = dz[photoz_inds[i]]/sqrt( Om*pow(1. + redshifts[i], 3) + (1. - Om)*pow(1. + redshifts[i], 3.*(1 + wDE)) );
                    }
                    if (cosmo_model == 5) {
                        dz_Hinv_term = dz[photoz_inds[i]]/sqrt( Om*pow(1. + redshifts_sort_fill[i], 3)
                           + (1. - Om)*pow(1. + redshifts_sort_fill[i], 3.*(1 + wDE + waDE))*exp(-3.*waDE*redshifts_sort_fill[i]/(1. + redshifts_sort_fill[i])) );
                    }
                }
                model_mu[i] = 5.*log10((1. + zhelio[i] + dz_term)*(r_com_sort[unsort_inds[i] + 1] + dz_Hinv_term)) + 43.22987755309658; //43.1586133146; h0=0.6774

                if (cosmo_model == 1) {
                    Hr[i] = sqrt( Om*pow(1. + redshifts[i], 3) + (1. - Om)) * r_com_sort[unsort_inds[i] + 1];
                }
                if (cosmo_model == 3) {
                    Hr[i] = sqrt( Om*pow(1. + redshifts_sort_fill[i], 3) + (1. - Om)*pow(1. + redshifts_sort_fill[i], 3.*(1 + wDE)) )* r_com_sort[unsort_inds[i] + 1];
                }
                if (cosmo_model == 5) {
                    Hr[i] = sqrt( Om*pow(1. + redshifts_sort_fill[i], 3)
                        + (1. - Om)*pow(1. + redshifts_sort_fill[i], 3.*(1 + wDE + waDE))*exp(-3.*waDE*redshifts_sort_fill[i]/(1. + redshifts_sort_fill[i])) )* r_com_sort[unsort_inds[i] + 1];
                }
            }
        }
        if (cosmo_model == 2) { // binned mu
            model_mu = dmu_dbin * mu_zbins + mu_const;
            for (i in 1:n_sne) {
                if (photoz_inds[i] > 0) {
                    // model_mu[i] = model_mu[i] + dz[photoz_inds[i]] * dmudz_dbin[i] * mu_zbins + mu_const;
                }
            }
        }

        if (cosmo_model == 6) { // binned comoving distance
            for (i in 1:n_zbins) {
                r_comove_bins[i] = 10^(0.2*(mu_zbins[i] - 43.1586133146))  /  (1. + zbins[i]);
            }

            model_mu = dmu_dbin * r_comove_bins;
            for (i in 1:n_sne) {
                model_mu[i] = 5.*log10((1. + zhelio[i])*model_mu[i]) + 43.1586133146;
            }
        }

        if (cosmo_model == 4) {
            for (i in 1:n_sne) {
                model_mu[i] = 5.*log10((1. + zhelio[i])*redshifts[i]/(1. + redshifts[i]) * (1. + (1./2.)*(1 - q0)*redshifts[i] - (1./6.)*(1. - q0 - 3.*q0*q0 + j0) * redshifts[i]*redshifts[i])
                   ) + 43.1586133146; // Equation 19 of Visser
            }
        }

        out[1]=model_mu;
        out[2]=Hr;
        // -------------End numerical integration---------------
        return out; 
    }

}


data {
    int<lower=0> n_sne; // number of SNe
    int<lower=0> n_samples;
    //int<lower=0> n_highz;
    int<lower=0> n_calib;
    int<lower=0> n_photoz; // number of SNe with photo-z's
    int n_x1c_star;    

    // int <lower=1, upper = n_samples> sample_list[n_sne];
    array[n_sne]int <lower=1, upper = n_samples> sample_list;
    array[n_sne]real <lower=0> redshifts;
    array[n_sne]real <lower=0> zhelio;
    matrix [n_sne, n_x1c_star] redshift_coeffs;
    

    int cosmo_model; // 1 => Om, 2 => Binned mu, 3 => Om-w, 4 => q0-j0, 5 => Om-w0-wa, 6 => comoving distance interpolation
    real fix_Om;
    int MB_by_sample;
    int <lower = 0, upper = 1> threeD_unexplained;

    int n_zbins;
    vector [n_zbins] zbins;
    vector [n_sne] mu_const;
    matrix [n_sne, n_zbins] dmu_dbin;
    matrix [n_sne, n_zbins] dmudz_dbin;
    vector [n_calib] z_low;
    // vector [n_sne] Hr;
    //vector [n_highz] z_high;


    array[n_sne] vector[3] obs_mBx1c;
    array[n_sne] matrix[3,3] obs_mBx1c_cov;
    array[n_sne] matrix[3, n_calib] d_mBx1c_d_calib;
    vector [n_sne] mass;
    vector [n_sne] mass_err;
    vector [n_sne] p_high_mass;

    int nzadd;
    array[2*(n_sne + nzadd) - 1]real redshifts_sort_fill ;
    array[n_sne + nzadd] int unsort_inds;

    int do_twoalphabeta;
    int do_host_mass;

    real outl_frac_prior_lnmean;
    real outl_frac_prior_lnwidth;

    vector [n_sne] mobs_cut0; // mobs_cut = mB + mobs_cut0 + mobs_cut1*c
    vector [n_sne] mobs_cut1;
    vector [n_samples] est_mobs_cuts;
    vector [n_samples] est_mobs_sigmas;

    array[n_photoz]vector [3] d_mBx1c_dz_list;
    vector [n_photoz] photo_z0;
    vector [n_photoz] photo_dz;
    vector [n_photoz] spike_redshift_prob; // E.g., 0.8
    vector [n_photoz] photo_spikez;
    array[n_sne] int <lower = 0, upper = n_sne> photoz_inds; // index of photo-z parameter (indexed from one) if photo-z, else 0
    //matrix [3,3] BAOCMB_Om_w0_wa_covmatrix;
    //vector [3] BAOCMB_Om_w0_wa_mean;
}

transformed data {
    int n_gauss = 4;
    vector [n_gauss] exp_approx_norm = [0.15038540936467037, 0.2993904768085472, 0.364279051173158, 0.18594506265362443]';
    vector [n_gauss] exp_approx_pos = [0.10329973984501734, 0.41080906196995237, 1.083137332416308, 2.427349566890827]';
    vector [n_gauss] exp_approx_width = [0.06596419371844692, 0.1910889454034621, 0.45516250820784515, 1.0637414822809306]';


    real a_ = 0.51;
    real K_ = 0.87;
    real deltaz_= 1/ (1/K_  +1);
    // vector [n_gauss] exp_approx_norm = [0.24410438, 0.43274856, 0.32314706]';
    // vector [n_gauss] exp_approx_pos = [0.16913558, 0.68695591, 1.9434773]';
    // vector [n_gauss] exp_approx_width = [0.11070724, 0.330062, 0.96505958]';

    print ("Version 1.71");

    // locs [ 0.16913558  0.68695591  1.9434773 ]
    // sigs [ 0.11070724  0.330062    0.96505958]
    // ampls [ 0.24410438  0.43274856  0.32314706]

    // locs [0.10329973984501734, 0.41080906196995237, 1.083137332416308, 2.427349566890827]
    // sigs [0.06596419371844692, 0.1910889454034621, 0.45516250820784515, 1.0637414822809306]
    // ampls [0.15038540936467037, 0.2993904768085472, 0.364279051173158, 0.18594506265362443]
}

parameters {
    vector [n_samples*MB_by_sample + 1*(1 - MB_by_sample)] MB;
    real <lower = -0.2, upper = 0.3> alpha_angle;
    real <lower = -1.4, upper = 1.4> beta_angle_red_low;

    real <lower = 0, upper = 1> Om;
    real <lower = 0> fs8_eff;
    real <lower=0> sigma_v;

    array[n_samples] real <lower=0.01, upper = 0.3> sigma_int;

    vector [n_sne] true_x1;
    vector [n_sne] true_cB;

    vector [n_calib] calibs_i;
}

transformed parameters {
    real alpha;
    real beta_B;
    vector [n_calib] calibs;

    calibs = calibs_i * fs8_eff; 
    alpha = tan(alpha_angle);
    beta_B = tan(beta_angle_red_low);
}

model {
    vector [3] Omw0wa_vect;

    // vector[n_sne] model_mu;
    array[2] vector[n_sne] model_mu_Hr; 
    array[n_sne] vector [3] model_mBx1c;
    array[n_sne] matrix [3,3] model_mBx1c_cov;

    array[n_samples] vector [3] sig_int_vector;
    array[3] vector [n_sne] sig_v;
    vector [n_sne] inl_loglike_by_SN;

    model_mu_Hr = calc_model_mu_Hr(n_sne, cosmo_model, nzadd, Om, redshifts, redshifts_sort_fill, unsort_inds, zhelio, photoz_inds);

    model_mBx1c_cov = obs_mBx1c_cov;

    // sigmaV part
    for (i in 1:n_sne) {
        //// NOTE THAT REALLY SHOULD BE DONE BY SAMPLE RATHER THAN REDSHIFT
        //// IT HAPPENS THAT THIS CONDITION BASED ON OBSERVED REDSHIFT DOES DISTINGUISH THE 2 SAMPLES BEING CONSIDERED
        if (redshifts[i] < 0.1) {
            sig_v[1][i]= (5/log(10.)) * (sigma_v / 299792.458) *(((1.+redshifts[i]) /model_mu_Hr[2][i]) - 1 );
        }
        else {
            sig_v[1][i]= (5/log(10.)) * (300 / 299792.458) *(((1.+redshifts[i]) /model_mu_Hr[2][i]) - 1 );
        }
        sig_v[2][i]= 0 ;
        sig_v[3][i]= 0 ;
    }

    for (i in 1:n_samples) {
        sig_int_vector[i][1] = sigma_int[i];        // This vector is in dispersion, not variance
        sig_int_vector[i][2] = 0.;
        sig_int_vector[i][3] = 0.;        
    }

    for (i in 1:n_sne) {
        for (j in 1:3) {
            model_mBx1c_cov[i][j,j] = model_mBx1c_cov[i][j,j] + sig_int_vector[sample_list[i]][j]^2 + sig_v[j][i]^2;
        }

        // model_mBx1c[i][1] = MB[1] + model_mu[i] - alpha*true_x1[i] + beta_B*true_cB[i]; 
        model_mBx1c[i][1] = MB[1] + model_mu_Hr[1][i] - alpha*true_x1[i] + beta_B*true_cB[i]; 
        model_mBx1c[i][2] = true_x1[i];
        model_mBx1c[i][3] = true_cB[i]; //+ true_cR[i];

        //// v TO m TRANSFORMATION DONE HERE
        inl_loglike_by_SN[i] = multi_normal_lpdf(obs_mBx1c[i] |
            model_mBx1c[i] + d_mBx1c_d_calib[i] * calibs * (5/log(10.)) / 299792.458 *(((1.+redshifts[i]) /model_mu_Hr[2][i]) - 1 ) , model_mBx1c_cov[i]);
    //                           sqrt(L) ev_i  * fs8_n z * dm/dv  
 
    }
                                        
    target += inl_loglike_by_SN;

    calibs_i ~ normal(0, 1);

    MB ~ normal(-19.12, 0.3);  

    //// PRIORS FOR DISPERSIONS
    // fs8_eff ~ cauchy(0,10);
    // sigma_v ~ cauchy(0,10);
    // sigma_int ~ cauchy(0,10);

    //// PRIORS FOR NON-FLAT SN PARAMETER DISTRIBUTIONS

    // target += log_sum_exp(log(a_-a_*deltaz_ + deltaz_) + normal_lpdf(true_x1 | 0.37, 0.61), log(1-deltaz_) + log(1-a_) + normal_lpdf(true_x1 | -1.22, 0.56));

}

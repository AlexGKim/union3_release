// Version History. Update in print statement in transformed data!
// Version 1.5. First official release with new selection-effect and population model!
// Version 1.6 (Sep-17-2020). Added full k-correction uncertainty propagation to selection-effect model.
// Version 1.7 (Jan-26-2024). Fixed parameter limits on selection-effect model outl_mBx1c_uncertainties.
// Version 1.71 (Sep-13-2024). Added lower limit to mobs_var_by_SN_except_c_R Thanks Aaron Do!


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
    vector [n_sne] Hr;
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
    // real <lower = -1.4, upper = 1.4> beta_angle_blue;
    real <lower = -1.4, upper = 1.4> beta_angle_red_low;
    // real <lower = -1.4, upper = 1.4> beta_angle_red_high;
    // real <lower = 9, upper = 11> mass_0;
    // real delta_0;
    // real <lower = 0, upper = 1> delta_h;

    real <lower = 0, upper = 1> Om;
    // real <lower = -2, upper = 0> wDE;
    // real <lower = -5, upper = 5> waDE;
    // real <lower = -2, upper = 2> q0;
    // real <lower = -5, upper = 5> j0;
    real <lower = 0> fs8_eff;
    real <lower=0> sigma_v;

    
    // vector [n_zbins] mu_zbins;

    array[n_samples] real <lower=0.01, upper = 0.3> sigma_int;
    // simplex [3] mBx1c_int_variance;

    vector [n_sne] true_x1;
    vector [n_sne] true_cB;
    //vector <lower = -0.25> [n_sne] true_cR_unit;


    // Population parameters:
  
    // vector <lower = -5, upper = 5> [n_x1c_star] x1_star;
    // vector <lower= 0.1, upper = 2> [n_x1c_star] R_x1;
    // vector <lower =-5, upper = 5> [n_x1c_star] tau_x1;

    // vector <lower = -0.5, upper = 0.5> [n_x1c_star] c_star;
    // vector <lower= 0.01, upper = 0.2> [n_x1c_star] R_c;
    // vector <lower = -0.5, upper = 0.5> [n_x1c_star] tau_c; // This doesn't actually have to be positive. If it's negative, the color distribution will have negative skew.

    // //

    vector [n_calib] calibs_i;

    

    // array[n_samples] real <lower = 100, upper = 200> mobs_cuts;
    // array[n_samples] real <lower = 0.1, upper = 3> mobs_cut_sigmas;

    // vector [n_photoz] dz;

    
    // real <lower = 0.2, upper = 2> outl_mBx1c_uncertainties_mB;
    // real <lower = 1, upper = 10> outl_mBx1c_uncertainties_x1;
    // real <lower = 0.1, upper = 1> outl_mBx1c_uncertainties_cB;
    //real <lower = 1, upper = 10> outl_mBx1c_uncertainties_cR_unit;
}

transformed parameters {
    array[n_sne] vector [3] model_mBx1c;
    array[n_sne] matrix [3,3] model_mBx1c_cov;
    // matrix [3,3] model_mBx1c_cov_outl [n_sne];

    array[n_samples] vector [3] sig_int_vector;
    // real <lower = -0.1, upper = 0.1> outl_frac;

    real alpha;
    real beta_B;
    // real beta_R_high;
    // real beta_R_low;
    //real <lower = 0. > sig_v;
    //real sig_v_m;
    array[3] vector [n_sne] sig_v;
    //vector [3,n_highz] sig_v[n_samples];

    real alpha_eff;
    real beta_eff;
    // real p_high_mass_eff;

    //vector [n_sne] true_cR;

    // vector [n_sne] x1_star_by_SN;
    // vector [n_sne] R_x1_by_SN;
    // vector [n_sne] tau_x1_by_SN;

    // vector [n_sne] c_star_by_SN;
    // vector [n_sne] R_c_by_SN;
    // vector [n_sne] tau_c_by_SN;
    // vector [n_calib] calibs_fs8;
    vector [n_calib] calibs;


    // vector [n_sne] mobs_by_SN_except_c_R;
    // vector <lower = 0.0001> [n_sne] mobs_var_by_SN_except_c_R; // Thanks Aaron Do!

    array[2*(n_sne + nzadd) - 1] real Hinv_sort_fill;
    array[n_sne + nzadd] real r_com_sort;
    vector [n_sne]  model_mu;

    // vector [n_sne] outl_loglike_by_SN;
    vector [n_sne] inl_loglike_by_SN;
    //vector [n_sne] PointPosteriors;

    // real this_MB;
    // real this_norm_LL;

    // vector [3] dz_deriv_term;
    // real dz_term;
    // real dz_Hinv_term;

    // vector [n_gauss] tmploglike_x1;
    //vector [n_gauss] tmploglike_c;

    // vector [n_zbins] r_comove_bins;

    // -------------Begin numerical integration-----------------


    

    if ((cosmo_model == 1) || (cosmo_model == 3) || (cosmo_model == 5)) {
        for (i in 1: 2*(n_sne + nzadd) - 1) {    // Inverse Hubble parameter
            if (cosmo_model == 1) {
                Hinv_sort_fill[i] = 1./sqrt( Om*pow(1. + redshifts_sort_fill[i], 3) + (1. - Om) );
            }
            // if (cosmo_model == 3) {
            //     Hinv_sort_fill[i] = 1./sqrt( Om*pow(1. + redshifts_sort_fill[i], 3) + (1. - Om)*pow(1. + redshifts_sort_fill[i], 3.*(1 + wDE)) );
            // }
            // if (cosmo_model == 5) {
            //     Hinv_sort_fill[i] = 1./sqrt( Om*pow(1. + redshifts_sort_fill[i], 3)
		    //                         + (1. - Om)*pow(1. + redshifts_sort_fill[i], 3.*(1 + wDE + waDE))*exp(-3.*waDE*redshifts_sort_fill[i]/(1. + redshifts_sort_fill[i])) );
            // }
        }

        // Integrate comoving r using Simpson's rule

        r_com_sort[1] = 0.; // Redshift = 0 should be first element!
        for (i in 2:(n_sne + nzadd)) {
            r_com_sort[i] = r_com_sort[i - 1] + (Hinv_sort_fill[2*i - 3] + 4.*Hinv_sort_fill[2*i - 2] + Hinv_sort_fill[2*i - 1])*(redshifts_sort_fill[2*i - 1] - redshifts_sort_fill[2*i - 3])/6.;
        }


        // for (i in 1:n_sne) {
        //     if (photoz_inds[i] == 0) {
        //         dz_term = 0.;
	    //     dz_Hinv_term = 0.;
        //     } else {
        //         dz_term = dz[photoz_inds[i]];
        //         if (cosmo_model == 1) {
        //             dz_Hinv_term = dz[photoz_inds[i]]/sqrt( Om*pow(1. + redshifts[i], 3) + (1. - Om) );
        //         }
        //         if (cosmo_model == 3) {
        //             dz_Hinv_term = dz[photoz_inds[i]]/sqrt( Om*pow(1. + redshifts[i], 3) + (1. - Om)*pow(1. + redshifts[i], 3.*(1 + wDE)) );
        //         }
        //         if (cosmo_model == 5) {
        //             dz_Hinv_term = dz[photoz_inds[i]]/sqrt( Om*pow(1. + redshifts_sort_fill[i], 3)
		//                                            + (1. - Om)*pow(1. + redshifts_sort_fill[i], 3.*(1 + wDE + waDE))*exp(-3.*waDE*redshifts_sort_fill[i]/(1. + redshifts_sort_fill[i])) );
        //         }
        //     }
        //     model_mu[i] = 5.*log10((1. + zhelio[i] + dz_term)*(r_com_sort[unsort_inds[i] + 1] + dz_Hinv_term)) + 43.1586133146;
        // }
        for (i in 1:n_sne) {
             //43.1586133146 ; # to get to H0=67.64
            model_mu[i] = 5.*log10((1. + zhelio[i] )*(r_com_sort[unsort_inds[i] + 1] )) + 43.229877553; 
        }

    }
    // if (cosmo_model == 2) { // binned mu
    //     model_mu = dmu_dbin * mu_zbins + mu_const;
    //     for (i in 1:n_sne) {
    //         if (photoz_inds[i] > 0) {
    //             // model_mu[i] = model_mu[i] + dz[photoz_inds[i]] * dmudz_dbin[i] * mu_zbins + mu_const;
    //         }
    //     }
    // }

    // if (cosmo_model == 6) { // binned comoving distance
    //     for (i in 1:n_zbins) {
    //         r_comove_bins[i] = 10^(0.2*(mu_zbins[i] - 43.1586133146))  /  (1. + zbins[i]);
    //     }
	
    //     model_mu = dmu_dbin * r_comove_bins;
	// for (i in 1:n_sne) {view.p
	//     model_mu[i] = 5.*log10((1. + zhelio[i])*model_mu[i]) + 43.1586133146;
    //     }
    // }

    // if (cosmo_model == 4) {
    //     for (i in 1:n_sne) {
    //         model_mu[i] = 5.*log10((1. + zhelio[i])*redshifts[i]/(1. + redshifts[i]) * (1. + (1./2.)*(1 - q0)*redshifts[i] - (1./6.)*(1. - q0 - 3.*q0*q0 + j0) * redshifts[i]*redshifts[i])
	//                            ) + 43.1586133146; // Equation 19 of Visser
    //     }
    // }


    // -------------End numerical integration---------------


    model_mBx1c_cov = obs_mBx1c_cov;
    // model_mBx1c_cov_outl = obs_mBx1c_cov;

    calibs = calibs_i * fs8_eff; 

    // sigmaV part
    for (i in 1:n_sne) {
        if (redshifts[i] < 0.1) {
            sig_v[1][i]= (5/log(10.)) * (sigma_v / 299792.458) *(((1.+redshifts[i]) * 299792.458/Hr[i]) - 1 ); }
        else {
             sig_v[1][i]= (5/log(10.)) * (300 / 299792.458) *(((1.+redshifts[i]) * 299792.458 /Hr[i]) - 1 );}
        sig_v[2][i]= 0 ;
        sig_v[3][i]= 0 ;}

    // outl_frac = 0;

    alpha = tan(alpha_angle);
    
    // if (do_twoalphabeta == 0) {
    //     beta_B = tan(beta_angle_red_low);
    //     beta_R_low = tan(beta_angle_red_low);
    //     beta_R_high = tan(beta_angle_red_low);
    // } else {
    //     beta_B = tan(beta_angle_blue);
    //     beta_R_low = tan(beta_angle_red_low);
    //     beta_R_high = tan(beta_angle_red_high);
    // }

    beta_B = tan(beta_angle_red_low);
    // beta_R_low = tan(beta_angle_red_low);
    // beta_R_high = tan(beta_angle_red_low);
    for (i in 1:n_samples) {
        // if (threeD_unexplained == 1) {
        //     sig_int_vector[i][1] = sqrt(mBx1c_int_variance[1])*sigma_int[i];        // This vector is in dispersion, not variance
        //     sig_int_vector[i][2] = sqrt(mBx1c_int_variance[2])*sigma_int[i]/0.14;
        //     sig_int_vector[i][3] = sqrt(mBx1c_int_variance[3])*sigma_int[i]/(-3.);
        // } else {
        //     sig_int_vector[i][1] = sigma_int[i];        // This vector is in dispersion, not variance
        //     sig_int_vector[i][2] = 0.;
        //     sig_int_vector[i][3] = 0.;
        // }
        sig_int_vector[i][1] = sigma_int[i];        // This vector is in dispersion, not variance
        sig_int_vector[i][2] = 0.;
        sig_int_vector[i][3] = 0.;        
        //mobs_cuts[i] = 21.5;
        //mobs_cut_sigmas[i] = 0.25;
    }



    for (i in 1:n_sne) {

        // if (do_host_mass == 1) {
    	//     p_high_mass_eff = (1.9*(1 - delta_h)/(1 + 0.9*exp(0.95*log(10.)*redshifts[i])) + delta_h)*p_high_mass[i];
        // } else {
        //     p_high_mass_eff = 0;
        // }

        for (j in 1:3) {
            model_mBx1c_cov[i][j,j] = model_mBx1c_cov[i][j,j] + sig_int_vector[sample_list[i]][j]^2 + sig_v[j][i]^2;
        }
	// model_mBx1c_cov_outl[i][1,1] = model_mBx1c_cov_outl[i][1,1] + outl_mBx1c_uncertainties_mB^2;


        // x1_star_by_SN[i] = dot_product(x1_star, redshift_coeffs[i]);
        // R_x1_by_SN[i] = dot_product(R_x1, redshift_coeffs[i]);
        // tau_x1_by_SN[i] = dot_product(tau_x1, redshift_coeffs[i]);

        // c_star_by_SN[i] = dot_product(c_star, redshift_coeffs[i]);
        // R_c_by_SN[i] = dot_product(R_c, redshift_coeffs[i]);
        // tau_c_by_SN[i] = dot_product(tau_c, redshift_coeffs[i]);

        //true_cR[i] = true_cR_unit[i]*tau_c_by_SN[i];


        // if (MB_by_sample == 1) {
        //     this_MB = MB[sample_list[i]];
        // } else {
        //     this_MB = MB[1];
        // }


	// mobs_by_SN_except_c_R[i] = this_MB + model_mu[i] + mobs_cut0[i] - alpha*x1_star_by_SN[i] + (beta_B + mobs_cut1[i])*c_star_by_SN[i] - delta_0*p_high_mass_eff;
	// mobs_var_by_SN_except_c_R[i] = mobs_cut_sigmas[sample_list[i]]^2
    //                                   + model_mBx1c_cov[i][1,1] + (mobs_cut1[i]*model_mBx1c_cov[i][3,3])^2 + 2.*mobs_cut1[i]*model_mBx1c_cov[i][1,3]
    //                                   + (alpha*R_x1_by_SN[i])^2 + ((beta_B + mobs_cut1[i])*R_c_by_SN[i])^2;


        model_mBx1c[i][1] = MB[1] + model_mu[i] - alpha*true_x1[i] + beta_B*true_cB[i]; //+(beta_R_low*(1 - p_high_mass_eff); 
                            //+ beta_R_high*p_high_mass_eff)*true_cR[i] - delta_0*p_high_mass_eff;
        model_mBx1c[i][2] = true_x1[i];
        model_mBx1c[i][3] = true_cB[i]; //+ true_cR[i];

	// model_mBx1c[i] = model_mBx1c[i]; //  + d_mBx1c_d_calib[i] * calibs;

    //     if (photoz_inds[i] == 0) {
    //         dz_deriv_term[1] = 0.;
    //         dz_deriv_term[2] = 0.;
    //         dz_deriv_term[3] = 0.;
    //     } else {
    //         dz_deriv_term = dz[photoz_inds[i]]*d_mBx1c_dz_list[photoz_inds[i]];
    //     }


        // for (g_ind in 1:n_gauss) { // gauss ind
        //     //tmploglike_c[g_ind]   = log(exp_approx_norm[g_ind]) + normal_log(true_cR_unit[i], exp_approx_pos[g_ind], exp_approx_width[g_ind]);
        //     tmploglike_x1[g_ind] = log(exp_approx_norm[g_ind]) + normal_log(true_x1[i], exp_approx_pos[g_ind]*tau_x1_by_SN[i] + x1_star_by_SN[i], sqrt((exp_approx_width[g_ind]*tau_x1_by_SN[i])^2 + R_x1_by_SN[i]^2));
        // }


	// outl_loglike_by_SN[i] = log(outl_frac) +
    //                     multi_normal_log(obs_mBx1c[i] + d_mBx1c_d_calib[i] * calibs + dz_deriv_term, model_mBx1c[i], model_mBx1c_cov_outl[i])
	// 				  + normal_log(true_x1[i], 0, outl_mBx1c_uncertainties_x1)
	// 				  + normal_log(true_cB[i], 0, outl_mBx1c_uncertainties_cB);
	// 				  //+ normal_log(true_cR_unit[i], 0, outl_mBx1c_uncertainties_cR_unit);


    //     this_norm_LL = 0.0001;
	// for (g_indx in 1:n_gauss) {
	//     for (g_indc in 1:n_gauss) {
	//         this_norm_LL += exp_approx_norm[g_indx]*exp_approx_norm[g_indc]*normal_cdf(   mobs_cuts[sample_list[i]], //  + d_mBx1c_d_calib[i][1] * calibs
	// 										      mobs_by_SN_except_c_R[i] + ((beta_R_low*(1 - p_high_mass_eff) + beta_R_high*p_high_mass_eff) + mobs_cut1[i])*exp_approx_pos[g_indc]*tau_c_by_SN[i] - alpha*exp_approx_pos[g_indx]*tau_x1_by_SN[i],
											   
    //                                                   	                                      sqrt(mobs_var_by_SN_except_c_R[i]
	// 										      + (((beta_R_low*(1 - p_high_mass_eff) + beta_R_high*p_high_mass_eff) + mobs_cut1[i])*exp_approx_width[g_indc]*tau_c_by_SN[i])^2
	// 										      + (alpha*exp_approx_width[g_indx]*tau_x1_by_SN[i])^2)   );
	//     }
	//

    
    // inl_loglike_by_SN[i] = log(1 - outl_frac)
    //                    + multi_normal_log(obs_mBx1c[i] + d_mBx1c_d_calib[i] * calibs + dz_deriv_term, model_mBx1c[i], model_mBx1c_cov[i])
	// 				  + normal_log(true_cB[i], c_star_by_SN[i], R_c_by_SN[i])
	// 				   + log_sum_exp(tmploglike_x1) // + log_sum_exp(tmploglike_c)

	//                                   + normal_cdf_log(mobs_cuts[sample_list[i]], //  + d_mBx1c_d_calib[i][1] * calibs
	// 				    obs_mBx1c[i][1] + d_mBx1c_d_calib[i][1] * calibs + mobs_cut0[i] + mobs_cut1[i]*(obs_mBx1c[i][3] + d_mBx1c_d_calib[i][3] * calibs),
	// 				    mobs_cut_sigmas[sample_list[i]])
    //                                       - log(this_norm_LL); //No calibration in this term, see above comment!
        inl_loglike_by_SN[i] = multi_normal_lpdf(obs_mBx1c[i] + d_mBx1c_d_calib[i] * calibs | model_mBx1c[i], model_mBx1c_cov[i]);
    // inl_loglike_by_SN[i] = multi_normal_lpdf(obs_mBx1c[i] + d_mBx1c_d_calib[i] * calibs + dz_deriv_term | model_mBx1c[i], model_mBx1c_cov[i]);
                      // + normal_log(true_cB[i], c_star_by_SN[i], R_c_by_SN[i])
                      //  + log_sum_exp(tmploglike_x1) // + log_sum_exp(tmploglike_c)

                      //                 + normal_cdf_log(mobs_cuts[sample_list[i]], //  + d_mBx1c_d_calib[i][1] * calibs
                      //   obs_mBx1c[i][1] + d_mBx1c_d_calib[i][1] * calibs + mobs_cut0[i] + mobs_cut1[i]*(obs_mBx1c[i][3] + d_mBx1c_d_calib[i][3] * calibs),
                      //   mobs_cut_sigmas[sample_list[i]])
                      //                     - log(this_norm_LL); //No calibration in this term, see above comment!
  }

}

model {
    vector [3] Omw0wa_vect;

    // ** Change by Alex : no outliers
    target += inl_loglike_by_SN;
    // for (i in 1:n_sne) {
    //     target += log_sum_exp(outl_loglike_by_SN[i], inl_loglike_by_SN[i]);
    // }

    // ** Change by Alex : no photoz
    // for (i in 1:n_photoz) {
    //     target += log_sum_exp(log(spike_redshift_prob[i]) + normal_log(dz[i], 0., 0.01),
    //                           log(1. - spike_redshift_prob[i]) + normal_log(dz[i], photo_z0[i] - photo_spikez[i], photo_dz[i]));
    // }

    calibs_i ~ normal(0, 1);

    //sig_v_i ~ normal(0, 1);
    
    // if (cosmo_model == 2) {
    //    MB ~ normal(-19, 0.001);
    // } else {
    

    // }

    // if (cosmo_model == 5) {
    //      Omw0wa_vect[1] = Om;   
    // 	 Omw0wa_vect[2] = wDE;
    // 	 Omw0wa_vect[3] = waDE;
	 
	//  //Omw0wa_vect ~ multi_normal(BAOCMB_Om_w0_wa_mean, BAOCMB_Om_w0_wa_covmatrix);
    // }

    MB ~ normal(-19, 0.3);
    // delta_0 ~ normal(0.0, 0.2);
    // mobs_cuts ~ normal(est_mobs_cuts, 0.5);
    // mobs_cut_sigmas ~ normal(est_mobs_sigmas, 0.25);


    if (fix_Om > 0) {
        Om ~ normal(fix_Om, 0.001);
    }

    // x1_star ~ normal(0., 2.);
    // tau_x1 ~ normal(-1., 2.);
    // R_x1 ~ normal(1, 2);

    // c_star  ~ normal(-0.1, 0.2);
    // tau_c ~ normal(0.1, 0.2);
    // R_c ~ normal(0.1, 0.2);
    // //true_cR ~ normal(0, 0.0001);
    

    // if (do_twoalphabeta == 0) {
    //    beta_angle_blue ~ normal(0, 1);
    // }

    //** Change by Alex commemnt out all outl stuff
    // outl_mBx1c_uncertainties_mB ~ normal(0.5, 0.5);
    // outl_mBx1c_uncertainties_x1 ~ normal(3, 3);
    // outl_mBx1c_uncertainties_cB ~ normal(0.5, 0.5);

    //outl_mBx1c_uncertainties_cR_unit ~ normal(10, 3);
    //sigma_int[1] ~ normal(0.02, 0.0001);

    //outl_frac ~ lognormal(outl_frac_prior_lnmean, outl_frac_prior_lnwidth);

    
    }

import multiprocessing
from matplotlib import use
# use("PDF")
import pickle
from numpy import *
import numpy as np

from cmdstanpy import CmdStanModel
import sys
import os
import time


from scipy.stats import scoreatpercentile
import helper_functions
from scipy.interpolate import interp1d
import gzip
from astropy.io import fits
from scipy.special import erf
import flip
import astropy.constants as acst
_C_LIGHT_KMS_ = acst.c.to('km/s').value
import astropy.cosmology as acosmo
import pandas as pd

def setup_multiprocessing_for_m1_m2():
    """Recommended setup for Mac M1/M2"""
    
    # Force spawn method for better compatibility
    if sys.platform == 'darwin':  # macOS
        multiprocessing.set_start_method('spawn')
        print("Using spawn method for Mac M1/M2 compatibility")
    
    # Get optimal process count
    cpu_count = multiprocessing.cpu_count()
    optimal_processes = min(cpu_count, 1)  # Don't exceed 8 for most tasks
    
    print(f"CPU cores: {cpu_count}")
    print(f"Recommended processes: {optimal_processes}")
    
    return optimal_processes


def main():
    optimal_processes = setup_multiprocessing_for_m1_m2()
    # optimal_processes = 1

    #######INPUT#########
    cosmo_model=1               
    stan_code_file = './stan_code_fs8_prune.stan' #your stan code file

    #number of iteration fro each chain and number of chains
    itera=500
    itera_warm = 500
    chains=4
    n_jobs = 4

    #txt file to save results and pickle file to save the chains
    # sys.stdout = open('res_fs8_snsim_lowz_svd_sqrt_vt_fitsigv_new_cov_errorpar.txt', 'w')
    fit_file='res_fs8_snsim_lowz_svd_sqrt_vt_fitsigv_new_cov_errorpar.pickle'

    #input cosmology for PV covariance
    cosmo_dic = {"h":0.6774, "omega_b":0.02230, "omega_cdm":0.1188,"sigma8":0.8159, "n_s":0.9667}#, 'mnu':0.0}
    cosmo= acosmo.FlatLambdaCDM(H0=67.74, Om0=0.3089)


    #input parameter for PV covariance
    sigma_u = 21
    kmin =  (2 * np.pi) / 2000 
    kmax = 0.5


    ########FUNCTIONS ###########

    #function to remove SN in the same host
    def find_samepos(df):
        samepos = []
        for ra in df.ra.unique():
            samera = df[df.ra == ra]
            for dec in samera.dec.unique():
                samedec = samera[samera.dec == dec]
                for zcos in samedec.zcos.unique():
                    samezcos = samedec[samedec.zcos == zcos]
                    if len(samezcos) >= 2:
                        samepos.append([df.index.get_loc(idx) for idx in samezcos.index.unique()])
        return samepos

    #function to compute PV power spectrum
    def init_PS(kmin, kmax,cosmo_dic):
      
        kh, _, _, ptt,fiducial = flip.power_spectra.compute_power_spectra(
        'class_engine',
        cosmo_dic, 
        0, kmin, 
        kmax, 
        1500, 
        normalization_power_spectrum="growth_rate",
        )
        return kh, ptt,fiducial


    #some functions copied from David code

    def add_zbins(stan_data, cosmo_model):
        # For binned mu
        
        stan_data["cosmo_model"] = cosmo_model

        print("min, max", stan_data["redshifts"].min(), stan_data["redshifts"].max())
        if stan_data["redshifts"].min() == stan_data["redshifts"].max() or ([2, 6].count(cosmo_model) == 0):
            stan_data["zbins"] = [stan_data["redshifts"][0]]
            stan_data["n_zbins"] = 1
            stan_data["dmu_dbin"] = ones([stan_data["n_sne"], stan_data["n_zbins"]], dtype=float64)
            stan_data["dmudz_dbin"] = zeros([stan_data["n_sne"], stan_data["n_zbins"]], dtype=float64)
            stan_data["mu_const"] = np.zeros(stan_data["n_sne"], dtype=np.float64)
            
            return stan_data


        assert cosmo_model == 6 or cosmo_model == 2
        zsort = np.sort(stan_data["redshifts"])

        print("zsort", zsort[-10:])

        zbins = [zsort[-1]*1.001]
        step = 10
        minstepsize = 0.1
        min_sn_bin = 10
        ind = -1 - min_sn_bin
        z_cutoff_for_05 = 0.8

        while step > minstepsize:
            step = zbins[0] - zsort[ind]
            minstepsize = ((zbins[0] + zsort[ind])*0.5 > z_cutoff_for_05)*0.05 + 0.05

            if step > minstepsize:
                zbins = [zsort[ind]] + zbins
                ind -= min_sn_bin

        print("zbins high z", zbins)


        zbins = np.concatenate((
            np.linspace(0.05, z_cutoff_for_05, int(np.around(z_cutoff_for_05/0.05))),
            np.linspace(z_cutoff_for_05, zbins[0], int(np.around((zbins[0] - z_cutoff_for_05)/0.1)) + 1)[1:-1],
            zbins))



        zbins = np.array(zbins)


        print("zbins", zbins, list(zbins))


        stan_data["zbins"] = zbins
        


    def get_redshifts(redshifts):
        appended_redshifts = arange(0., 2.51, 0.1)
        tmp_redshifts = concatenate((redshifts, appended_redshifts))
        
        sort_inds = list(argsort(tmp_redshifts))
        unsort_inds = [sort_inds.index(i) for i in range(len(tmp_redshifts))]
        
        tmp_redshifts = sort(tmp_redshifts)
        redshifts_sort_fill = sort(concatenate((tmp_redshifts, 0.5*(tmp_redshifts[1:] + tmp_redshifts[:-1]))))
        
        return redshifts_sort_fill, unsort_inds, len(appended_redshifts)


    def get_redshift_coeffs(z_list, p_high_mass, separate_mass_x1c, redshift_coeff_type,sample_list):
        """redshift_coeff_type could be ("a", 1) or ("a", 3) for a population that varies with a(t)
        redshift_coeff_type could be ("sample", [0.0, 0.4, 1.0]) for a population that is allowed to be different low-z, mid-z, high-z"""
        
        if redshift_coeff_type[1].count("."):
            n_z = len(redshift_coeff_type[1:])
        else:
            n_z = int(redshift_coeff_type[1])
        
        actual_n_x1c_star = n_z*(1 + separate_mass_x1c)

        redshift_coeffs = np.zeros([len(z_list), actual_n_x1c_star], dtype=np.float64)

        if n_z == 1:
            if separate_mass_x1c:
                redshift_coeffs[:,0] = p_high_mass
                redshift_coeffs[:,1] = 1 - p_high_mass
            else:
                redshift_coeffs += 1

            return redshift_coeffs

        if redshift_coeff_type[0] == "a":
            a_list = 1./(1. + np.array(z_list))
            a_nodes = np.linspace(min(a_list) - 1e-5, max(a_list) + 1e-5, n_z)
        
            for i in range(len(z_list)):
                for j in range(n_z):
                    coeffs = zeros(n_z, dtype=float64)
                    coeffs[j] = 1

                    ifn = interp1d(a_nodes, coeffs, kind = 'linear')

                    if separate_mass_x1c:
                        redshift_coeffs[i,j] = ifn(a_list[i])*p_high_mass[i]
                        redshift_coeffs[i,n_z + j] = ifn(a_list[i])*(1. - p_high_mass[i])
                    else:
                        redshift_coeffs[i,j] = ifn(a_list[i])


        elif redshift_coeff_type[0] == "sample":
            zs_to_match = np.array([float(item) for item in redshift_coeff_type[1:]])
            print("zs_to_match", zs_to_match)
            
            for set_ind in np.unique(set_list):
                mean_z = np.mean(z_list[np.where(set_list == set_ind)])
                
                j = np.argmin(np.abs(zs_to_match - mean_z))
                print("set_ind", set_ind, "mean_z", mean_z, "j", j)

                if separate_mass_x1c:
                    redshift_coeffs[:,j] += (set_list == set_ind)*p_high_mass
                    redshift_coeffs[:,n_z + j] += (set_list == set_ind)*(1. - p_high_mass)
                else:
                    redshift_coeffs[:,j] += (set_list == set_ind)*1.
        else:
            assert 0, "Unknown redshift_coeff_type " + str(redshift_coeff_type)
        
        return redshift_coeffs




    ######## MAIN #######

    ####read data and prepare stan_data#########

    ddf=pd.read_parquet('snsim_highz.parquet')
    df=pd.read_parquet('snsim_lowz.parquet')


    #remove SN in the same host for the lowz sample
    samepos_list = find_samepos(df)
    print(len(samepos_list))

    same_pos_mask = np.ones(len(df), dtype='bool')
    for same_pos in samepos_list:
        same_pos_mask[same_pos[1:]] = False

    wfd = df[same_pos_mask]


    #prepare the data for stan
    NSN =len(np.append(wfd.zobs.values,ddf.zobs.values))
    sample_list=np.append(np.asarray([1 for i in range(len(wfd))]),np.asarray([2 for i in range(len(ddf))]))
    zz=np.append(wfd.zobs.values,ddf.zobs.values)

    redshifts_sort_fill, unsort_inds, nzadd = get_redshifts(zz)

    mass=np.append(np.zeros(len(wfd))+11,np.zeros(len(ddf))+10)
    mass_err=np.ones(NSN)
    p_high_mass = 0.5*(1. + erf((mass - 10.)/(np.sqrt(2.) * mass_err)))
    redshift_coeffs= get_redshift_coeffs(zz, np.ones(NSN), 0, ['a', '1'],sample_list)


    #mb-x1-c data + covariance

    emb = np.random.normal(loc=0.0, scale=1.e-4, size=len(wfd))
    ec = np.random.normal(loc=0.0, scale=1.e-4, size=len(wfd))
    ex1 = np.random.normal(loc=0.0, scale=1.e-4, size=len(wfd))

    obs_mBx1c_wfd=[]
    for i in range(len(wfd)):
        obs_mBx1c_wfd.append([wfd.mb.values[i] + emb[i] , wfd.x1.values[i]+ex1[i], wfd.c.values[i]+ec[i]])
    obs_mBx1c_ddf=[]


    emb1 =np.random.normal(loc=0.0, scale=1.e-4, size=len(ddf))
    ec1 =np.random.normal(loc=0.0, scale=1.e-4, size=len(ddf))
    ex11 = np.random.normal(loc=0.0, scale=1.e-4, size=len(ddf)) 

    for i in range(len(ddf)):
        obs_mBx1c_ddf.append([ddf.mb.values[i]+ emb1[i], ddf.x1.values[i]+ ex11[i], ddf.c.values[i]+ec1[i]])
        
    obs_mBx1c=np.concatenate([obs_mBx1c_wfd,obs_mBx1c_ddf],axis=0)


    #####build covarances#########


    mBx1c_cov=np.zeros([NSN,3,3], dtype=float64)


    for i in range(NSN):
        mBx1c_cov[i][0][0] = 1.e-8
        mBx1c_cov[i][1][1] = 1.e-8
        mBx1c_cov[i][2][2] = 1.e-8


        

    #flip covariance###
    kh, ptt,fiducial = init_PS(1e-5, kmax,cosmo_dic)
    fs8_fid = fiducial["fsigma_8"]


    pw_dic = {'vv': [[kh, ptt * flip.utils.Du(kh, sigma_u)**2]]}
    rcom = cosmo.comoving_distance(wfd.zobs.values).value * cosmo.h

    COV = flip.covariance.CovMatrix.init_from_flip('carreres23', "velocity", pw_dic, 
                                                    coordinates_velocity=[wfd.ra.values, wfd.dec.values, rcom], 
                                                    kmin=kmin, number_worker=optimal_processes,)



    # #from CVV to Cmm

    # Hr = cosmo.H(wfd.zobs.values).value * cosmo.comoving_distance(wfd.zobs.values).value
    # vcoeff = np.diag(-5 / (np.log(10) * _C_LIGHT_KMS_) * ((1 + wfd.zobs.values) * _C_LIGHT_KMS_ / Hr - 1))


    vel_cov = COV.compute_covariance_sum({'fs8':1,'sigv':0},np.zeros(len(wfd)))
    # COV_mm = np.array(vcoeff.T @ vel_cov @ vcoeff)
    # _,s,v=np.linalg.svd(COV_mm)
    # vt = v.T
    #eigenvactors of pv covariance

    u,s,vt=np.linalg.svd(vel_cov,hermitian=True)
    d_mBx1c_dcalib=np.zeros([NSN,3,len(wfd)], dtype=float64)
    dum=np.zeros([NSN,3,len(wfd)], dtype=float64)
    for i in range(len(wfd)):
        d_mBx1c_dcalib[i][:len(wfd)] = u[i]*np.sqrt(s)

    # for i in range(len(wfd)):
    #     val = vt[i] * np.sqrt(s[i])
    #     for j in range(len(val)):
    #         dum[i][0][j] = val[j]
    

    ####### stan_data#########

    stan_data =  {   "n_sne": NSN, 
                     "nzadd": nzadd,
                     "n_samples": 2,
                     "n_highz" : len(ddf),
                     "redshift_coeffs": redshift_coeffs,
                     "z_low": wfd.zobs.values,
                     "z_high": ddf.zobs.values,
                     # "Hr" : cosmo.H(zz).value * cosmo.comoving_distance(zz).value,
                     "n_x1c_star": len(redshift_coeffs[0]), # 3 = 3 scale-factor nodes
                     "threeD_unexplained": 0,
                     "mass": mass,#+10 ,
                     "mass_err": mass_err,
                     "p_high_mass": p_high_mass,
                     "do_host_mass": 0, 
                     "fix_Om": 0, 
                     "MB_by_sample": 0, 
                     # The +1 here is for Stan's indexing, which is from 1 not 0
                     "sample_list": sample_list,
                     "zhelio": zz,
                     "redshifts":zz,
                     "redshifts_sort_fill": redshifts_sort_fill,
                     "unsort_inds": unsort_inds,
                     
                     "do_blind": 0,
                     "do_twoalphabeta": 0 ,

                     "outl_frac_prior_lnmean": 0.0001, #as in Union3
                     "outl_frac_prior_lnwidth": 0.0001,

                     "n_calib": len(wfd), #here the velocities
                     "d_mBx1c_d_calib": d_mBx1c_dcalib,
                     "d_mBx1c_dz_list":  np.zeros([0,3], dtype=float64),
                     "obs_mBx1c": obs_mBx1c,
                     "obs_mBx1c_cov": mBx1c_cov,
                  
                     "n_photoz": 0,
                     "photo_z0": [],
                     "photo_dz": [],
                     "spike_redshift_prob": [0.8]*0,
                     "photoz_inds": np.zeros(NSN,dtype=int8),
                     "photo_spikez": [],


                     "est_mobs_cuts": [150,150],
                     "est_mobs_sigmas": [1,1],
                     "mobs_cut0": np.zeros(NSN), 
                     "mobs_cut1": np.zeros(NSN)
                     #"BAOCMB_Om_w0_wa_mean": BAOCMB_Om_w0_wa_mean, "BAOCMB_Om_w0_wa_covmatrix": BAOCMB_Om_w0_wa_covmatrix
                 }

    stan_data = add_zbins(stan_data, cosmo_model)

    print("nzadd ", stan_data['nzadd'])


    ################################################# Init FN ###################################################

    #init of the chains

    def init_fn():
        n_sne = stan_data["n_sne"]
        n_samples = stan_data["n_samples"]
        print("n_sne ", n_sne)
        print("n_samples ", n_samples)

        if stan_data["cosmo_model"] == 2 or stan_data["cosmo_model"] == 6:
            zbins_tmp = np.array(stan_data["zbins"])
            mu_init = 43.2 + 5*np.log10((zbins_tmp - 0.225*zbins_tmp**2.)*(1. + zbins_tmp))
        
        else:
            mu_init = np.zeros(stan_data["n_zbins"], dtype=np.float64)
            
                
        return {"MB": np.random.random(size = [(n_samples - 1)*stan_data["MB_by_sample"] + 1])*0.2 - 19.1,
                "Om": 0.3,
                "wDE": -1.01,
                "fs8_eff": 0.95,
                "mu_zbins": mu_init,
                "alpha_angle": np.arctan(random.random()*0.2),
                "beta_angle_blue": np.arctan(random.random()*0.5 + 2.5),
                "beta_angle_red_low": np.arctan(random.random()*0.5 + 2.5),
                "beta_angle_red_high": np.arctan(random.random()*0.5 + 2.5),
                #"log10_sigma_int": log10(random.random(size = n_samples)*0.1 + 0.1),
                "mBx1c_int_variance": [0.9, 0.05, 0.05],
                #"mass_0": 10,
                "delta_0": np.random.random()*0.05,
                "delta_h": 0.5,
                "calibs": np.random.normal(size = len(wfd))*0.01,
                #"blind_values": [0.]*n_samples,
                
                "true_cB": np.random.random(size = n_sne)*0.02 - 0.01 + clip(np.append(wfd.c.values,ddf.c.values)/2., -0.2, 1.0),
                "true_cR_unit": np.random.random(size = n_sne)*0.5 + 0.5, #random.random(size = n_sne)*0.01 + clip(the_data["c_list"]/2., 0, 1.0),
                "true_x1": np.random.random(size = n_sne)*0.2 - 0.1 + np.append(wfd.x1.values,ddf.x1.values),

                "x1_star": np.random.random(size = stan_data["n_x1c_star"])*0.5,
                "tau_x1": -np.random.random(size = stan_data["n_x1c_star"]),
                "R_x1": np.random.random(size = stan_data["n_x1c_star"])*0.5 + 0.25,

                "c_star": -np.random.random(size = stan_data["n_x1c_star"])*0.05,
                "tau_c": np.random.random(size = stan_data["n_x1c_star"])*0.05,
                "R_c": np.random.random(size = stan_data["n_x1c_star"])*0.05 + 0.02,
                #"sigma_int":[0.03,0.06],
                "outl_frac": np.random.random()*0.02 + 0.01,
                "mobs_cuts": stan_data["est_mobs_cuts"] + np.random.normal(size = n_samples)*0.1, 
                "mobs_cut_sigmas": [0.5]*n_samples,
                

                "dz": np.random.normal(size = stan_data["n_photoz"])*0.01
            }
                
                
    #####RUN STA N####
        
    print("Running...")

    # Save timestamp
    start = time.time()



    # sm = pystan.StanModel(file=stan_code_file)
    # fit = sm.sampling(data=stan_data,
    #                   iter=itera, chains=chains, n_jobs = n_jobs, refresh = 20, 
    #                   init = init_fn, control=dict(max_treedepth=11,adapt_delta = 0.85))

    model = CmdStanModel(stan_file=stan_code_file)
    # fit = model.sample(data=stan_data,
    #                   iter=itera, chains=chains, n_jobs = n_jobs, refresh = 20, 
    #                   init = init_fn, control=dict(max_treedepth=11,adapt_delta = 0.85))

    fit = model.sample(data=stan_data,
             iter_sampling=itera,
             iter_warmup=itera_warm,
             chains=chains,
             parallel_chains=n_jobs,
             # refresh=20,
             inits=init_fn(),
             max_treedepth=11,
             adapt_delta=0.85, show_progress=True)

    df = fit.draws_pd()
    df.to_pickle("/Users/akim/Projects/union3_release/output/result.pkl")

    fit.save_csvfiles(dir="/Users/akim/Projects/union3_release/output")

    # #save results 
    # fit_params = fit.extract(permuted = True)
    # pickle.dump(fit_params, gzip.open(fit_file, "wb"))
    # try:
    #     fit_params = filter_fit_params(fit_params, "MB", chains, itera/2) # burns the first half of the chain, so iter/2
    # except:
    #     print("Couldn't filter bad chains! One or more chains may be bad!")



    # try:
    #     mu_cov = np.cov(fit_params["mu_zbins"].T)
    #     whole_mat = np.zeros([len(mu_cov) + 1]*2, dtype=np.float64)
    #     whole_mat[1:, 1:] = np.linalg.inv(mu_cov)
    #     whole_mat[1:, 0] = np.median(fit_params["mu_zbins"], axis = 0)
    #     whole_mat[0, 1:] = stan_data["zbins"]
    #     np.save('mu_cov.npy',whole_mat)
         
    # except:
    #     print("Couldn't save whole_mat")


    # del_keys = []
    # for key in fit_params:
    #     sh = np.array(fit_params[key].shape)

    #     if np.any(sh[1:] > 10000):
    #         print(key, " is too big to save!", sh)
    #         del_keys.append(key)

    # print("del_keys", del_keys)
    # for key in del_keys:
    #     del fit_params[key]

        
    # pickle.dump(fit_params, gzip.open(fit_file, "wb"))
    # print("I hope you have a log file:")

    try:
        print(fit.summary())
    except:
        print("Couldn't print fit! Something is very wrong!")


    # Save timestamp
    end = time.time()
    print('running time')
    print(end - start)

if __name__ == "__main__":
    main()
"""
Adjusted from ep-stan package
"""

from __future__ import division
import os, argparse, time
import numpy as np
from util import load_stan, suppress_stdout
import pystan

CUR_PATH = os.path.dirname(os.path.abspath(__file__))
MOD_PATH = CUR_PATH
RES_PATH = os.path.join(CUR_PATH, 'results')
# Temp fix for the RandomState seed problem with pystan in 32bit Python.
# Detect automatically if in 32bit mode
TMP_FIX_32BIT = os.sys.maxsize <= 2**32

MC_FULL_OPT = dict(
        chains  = 4,
        iter    = 1000,
        warmup  = 500,
        thin    = 2,
    )

def gmm_sampling(X, y, J, prior_mu, prior_precision, w, std_noise, seed_mcmc=0, \
        mc_full_opt=None, save_res=False, ids=None):
    """
    Sample the means of the components from the posterior
    """
    model_name = 'gmm2'
    #print "Full model {} ...".format(model_name)
        
    seed = np.random.RandomState(seed=seed_mcmc)
    # Temp fix for the RandomState seed problem with pystan in 32bit Python
    seed = seed.randint(2**31-1) if TMP_FIX_32BIT else seed
        
    # compute prior_cov using prior_precision  
    prior_cov = []
    if len(prior_precision.shape) == 2:	# diagonal matrix
        for j in xrange(J):
            prior_cov.append(np.diag(1.0 / prior_precision[j]))
    else:
        for j in xrange(J):
            prior_cov.append(np.linalg.inv(prior_precision[j]))
    prior_cov = np.array(prior_cov)
    
    # sampler settings
    if mc_full_opt is None:
        mc_full_opt = MC_FULL_OPT

    # construct data for stan model
    data = dict(        
        D = X.shape[1],
        N = X.shape[0],
        J = J,
        X = X,
        weights = w, 
        mu = prior_mu,
        sigma = prior_cov,
        std_noise = std_noise
    )
    if model_name == 'gmm2':
        data['y'] = np.array(y+1, dtype = int)
    
    # load stan model
    stan_model = load_stan(os.path.join(MOD_PATH, model_name))
        
    # sample and extract parameters
    with suppress_stdout():
        fit = stan_model.sampling(
            data = data,
            seed = seed,
            **mc_full_opt
        )
    samp = fit.extract(pars='phi')['phi']
    m_phi_full = samp.mean(axis=0)
    var_phi_full = []
    for j in xrange(J):
        var_phi_full.append(np.cov(samp[:, j], rowvar = 0))
    var_phi_full = np.array(var_phi_full)
    samp_phi = np.copy(samp)
        
    # Save results
    if save_res:
        if not os.path.exists(RES_PATH):
            os.makedirs(RES_PATH)
        if ids:
            filename = 'res_f_{}_{}.npz'.format(model_name, ids)
        else:
            filename = 'res_f_{}.npz'.format(model_name)
        np.savez(
            os.path.join(RES_PATH, filename),
            m_phi_full   = m_phi_full,
            var_phi_full = var_phi_full,
            samp_phi_full= samp_phi
        )
        print "Full model results saved."
    
    return m_phi_full, var_phi_full, samp_phi


"""A simulated experiment model used by the sckript fit.py

Model name: m1b
Definition:
    group index j = 1 ... J
    input index d = 1 ... D
    explanatory variable x = [x_1 ... x_D]
    response variable y
    shared parameter beta = [[beta_11 ... beta_1D] ... [beta_J1 ... beta_JD]]
    local parameter x = [x_1, ..., x_N]
    x ~ multinomial(w_1, ..., w_J)
    y ~ normal(beta_x, sigma_yH * I)
    beta_j ~ N(mu_bjH, sigma_bjH), for all d
    phi = [beta]

"""

# Licensed under the 3-clause BSD license.
# http://opensource.org/licenses/BSD-3-Clause
#
# Copyright (C) 2014 Tuomas Sivula
# All rights reserved.

from __future__ import division
import numpy as np
from common import data, calc_input_param_classification
from scipy.stats import norm
import pickle


# ------------------------------------------------------------------------------
# >>>>>>>>>>>>> Configurations start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ------------------------------------------------------------------------------

# ====== Model parameters ======================================================
# If BETA is None, it is sampled from N(0,SIGMA_B)
BETA = None

# ====== Prior =================================================================
# Prior for beta
M0_B = 0
V0_B = 1.5**2
P_X = 1	# TODO: need to normalise it later

# ====== Regulation ============================================================
# Min for abs(sum(beta))
B_ABS_MIN_SUM = 1e-4

# ------------------------------------------------------------------------------
# <<<<<<<<<<<<< Configurations end <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ------------------------------------------------------------------------------


class model(object):
    """Model definition.
    
    Parameters
    ----------
    J : int
        Number of groups
    
    D : int
        Number of inputs
    
    npg : {int, seq of ints}
        Number of observations per group (constant or [min, max])
    
    """
    
    def __init__(self, J, D, sigma_noise, MU_B = 0, SIGMA_B = 1.0, w = 1.0):
        self.J = J
        self.D = D
        self.dphi = D
        self.sigma_noise = sigma_noise
        self.MU_B = MU_B
        self.SIGMA_B = SIGMA_B
        if type(w) == float:
            self.w = np.ones(self.J) / float(self.J)
        else:
            self.w = np.array(w)
            self.w /= float(self.w.sum())
    
    def simulate_data(self, num_data, seed=None):
        """Simulate data from the model.
        
        Returns models.common.data instance
        
        """
        # Localise params
        J = self.J
        D = self.D
        MU_B = self.MU_B
        SIGMA_B = self.SIGMA_B
        N = num_data
        
        # Set seed
        rnd_data = np.random.RandomState(seed=seed)
        
        if BETA is None:
            beta = MU_B + rnd_data.randn(J, D)*SIGMA_B
        else:
            beta = BETA
        
        # Regulate beta
        for j in xrange(J):
            beta_sum = np.sum(beta[j])
            while np.abs(beta_sum) < B_ABS_MIN_SUM:
                # Replace one random element in beta
                index = rnd_data.randint(D)
                beta_sum -= beta[j,index]
                beta[j,index] = MU_B + rnd_data.randn()*SIGMA_B	# TODO
                beta_sum += beta[j,index]
        
        phi_true = beta
        
        X = np.empty(N, dtype=int)
        for n in xrange(N):
            X[n] = rnd_data.choice(J, p = self.w)	# assign to the jth group
        y = beta[X] + rnd_data.randn(N, D)*self.sigma_noise
        
        tmp = X
        X = y
        y = tmp
        
        data_gen = {'X': X, 'y': y, 'mu_beta': MU_B, 'sigma_beta': SIGMA_B, \
            'sigma_noise': self.sigma_noise, 'beta': beta}
        
        return data_gen
    
    def load_data(self, name, seed=None):
        """Simulate data from the model.
        
        Returns models.common.data instance
        
        """
        
        # load data
        data_full = np.loadtxt('data/%s.txt' % name)

        # We obtain the features and the targets

        X = data_full[ :, range(data_full.shape[ 1 ] - 1) ]
        y = data_full[ :, data_full.shape[ 1 ] - 1 ].astype(int)

        # We create the train and test sets with 90% and 10% of the data

        permutation = np.random.choice(range(X.shape[ 0 ]),
            X.shape[ 0 ], replace = False)
        size_train = np.round(X.shape[ 0 ] * 0.9)
        index_train = permutation[ 0 : size_train ]
        index_test = permutation[ size_train : ]

        X_train = X[ index_train, : ]
        y_train = y[ index_train ]
        X_test = X[ index_test, : ]
        y_test = y[ index_test ]
        
        # Localise params
        J = self.J
        self.D = X_train.shape[1]
        D = self.D
        self.npg = int(X_train.shape[0] / float(J))
        npg = self.npg
        
        X_train = X_train[:(self.npg*self.J), :]
        y_train = y_train[:(self.npg*self.J)]
        
        # Set seed
        rnd_data = np.random.RandomState(seed=seed)
        
        # Parameters
        # Number of observations for each group
        if hasattr(npg, '__getitem__') and len(npg) == 2:
            Nj = rnd_data.randint(npg[0],npg[1]+1, size=J)
        else:
            Nj = npg*np.ones(J, dtype=np.int64)
        # Total number of observations
        N = np.sum(Nj)
        # Observation index limits for J groups
        j_lim = np.concatenate(([0], np.cumsum(Nj)))
        # Group indices for each sample
        j_ind = np.empty(N, dtype=np.int64)
        for j in xrange(J):
            j_ind[j_lim[j]:j_lim[j+1]] = j
        
        if BETA is None:
            beta = rnd_data.randn(D)*SIGMA_B
        else:
            beta = BETA
        
        # Regulate beta
        beta_sum = np.sum(beta)
        while np.abs(beta_sum) < B_ABS_MIN_SUM:
            # Replace one random element in beta
            index = rnd_data.randint(D)
            beta_sum -= beta[index]
            beta[index] = rnd_data.randn()*SIGMA_B
            beta_sum += beta[index]
        
        phi_true = beta
        # Determine suitable mu_x and sigma_x
        mu_x_j, sigma_x_j = calc_input_param_classification(np.zeros(J), beta)
        
        return data(
            X_train, y_train, {'mu_x':mu_x_j, 'sigma_x':sigma_x_j}, y_train, Nj, j_lim, 
            j_ind, {'phi':phi_true, 'beta':beta}
        ), X_test, y_test
    
    def get_prior(self):
        """Get prior for the model.
        
        Returns: S, m, Q, r
        
        """
        D = self.D
        # Moment parameters of the prior (transposed in order to get
        # F-contiguous)
        S0 = np.diag(np.ones(D)*V0_B).T
        m0 = np.ones(D)*M0_B
        # Natural parameters of the prior
        Q0 = np.diag(np.ones(D)/V0_B).T
        r0 = np.ones(D)*(M0_B/V0_B)
        p_X = np.ones(self.J) * P_X
        p_X /= float(p_X.sum())
        
        return S0, m0, Q0, r0#, p_X
    
    def get_param_definitions(self):
        """Return the definition of the inferred parameters.
        
        Returns
        -------
        names : seq of str
            Names of the parameters
        
        shapes : seq of tuples
            Shapes of the parameters
        
        hiers : seq of int 
            The indexes of the hierarchical dimension of the parameter or None
            if it does not have one.
        
        """
        names = ['beta']
        shapes = [self.D]
        hiers = [None]
        return names, shapes, hiers



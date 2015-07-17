import numpy as np
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy.misc import logsumexp

def comp_update(mean, cov, r, MU, SIG, sig1, average = 'none', return_mean = False):
	# mean = r * (y - MU)
	# cov = E[ r *(y - MU) ** 2] - mean ** 2
	tmp = SIG / (SIG + sig1)
	mu = MU + tmp * mean
	sig = SIG + tmp ** 2 * cov - r * tmp * SIG
	if return_mean:
		return mu, sig
	if average == 'moments':
		#second_moment = sig + mu ** 2
		#mu = mu.mean(0)
		#second_moment = second_moment.mean(0)
		#sig = second_moment - mu ** 2
		mu = mu.mean(0)
		sig = sig.mean(0)
	inv_sig = 1 / sig - 1 / SIG
	mu_inv_sig = mu / sig - MU / SIG
	if average == 'params':
		mu_inv_sig = mu_inv_sig.mean(0); inv_sig = inv_sig.mean(0)
	return mu_inv_sig, inv_sig
	
def gmm_updates(data, Var, mean, approx_x = False, pred = False, full_cov = False):
    """
    Analytical solution for EP on the means of Gaussian mixtures.
    Var, mean is the paramter of mean cavity, w contains prior for 
    the components of mixture. noise_sig is the variance of noise
    dim(Var) = [J, D, D], dim(mean) = [J, D]
    data = (y, w, noise_sig)
    dim(w) = J
    """
    # unpack data inputs
    y = data[0]; w = data[1]; noise_sig = data[2]
    
    if y.shape[-1] == 1:
        Var = Var[:, 0]
        mean = mean[:, 0]
    
    # compute r_k
    K = w.shape[0]
    D = y.shape[0]
    precision_tilted = []
    r_k = []
    for k in xrange(K):
        if full_cov is True:
            #print full_cov
            Var_tilted = Var[k] + np.eye(D) * noise_sig	# matrix version
            precision_tilted.append(np.linalg.inv(Var_tilted))	# matrix version
            #pdf = multivariate_normal.pdf(y, mean = mean[k], cov = Var_tilted)	# matrix version
            m_diff = y - mean[k]
            log_pdf = -0.5 * np.dot(m_diff, np.dot(precision_tilted[-1], m_diff))
            log_pdf -= 0.5 * (np.log(2 * np.pi) * D + np.log(np.linalg.det(Var_tilted)))
        else:
            #print full_cov
            Var_tilted = Var[k] + noise_sig	# diagonal matrix
            precision_tilted.append(1.0 / Var_tilted)	# diagnoal matrix version
            #pdf = norm.pdf(y, loc = mean[k], scale = Var_tilted).prod()
            m_diff = y - mean[k]
            log_pdf = -0.5 * m_diff ** 2 / Var_tilted
            log_pdf -= 0.5 * (np.log(2 * np.pi) + np.log(Var_tilted))
            log_pdf = log_pdf.sum()
        # take log
        r_k.append(log_pdf + np.log(w[k]))        

    r_k = np.array(r_k)
    logZ = logsumexp(r_k)
    r_k = np.exp(r_k - logZ)
    if pred == True:
        return r_k, logZ
    
    precision_tilted = np.array(precision_tilted)

    # compute updates (forward)
    MU = []; SIG = []
    for k in xrange(K):
        if full_cov:
            dot_tmp = np.dot(Var[k], precision_tilted[k])	# matrix version
            dot_tmp2 = np.dot(dot_tmp, (y - mean[k]))	# matrix versioin
            s = Var[k] + r_k[k] * (1 - r_k[k]) * np.outer(dot_tmp2, dot_tmp2)	# matrix version
            s = s - r_k[k] * np.dot(dot_tmp, Var[k])	# matrix version
        else:
            dot_tmp = Var[k] * precision_tilted[k]	# diagonal matrix version
            dot_tmp2 = dot_tmp * (y - mean[k])	# diagonal matrix version
            s = Var[k] + r_k[k] * (1 - r_k[k]) * (dot_tmp2 ** 2)	# diagonal matrix version
            s = s - r_k[k] * (dot_tmp * Var[k])	# diagonal matrix version 
        m = mean[k] + r_k[k] * dot_tmp2    
        MU.append(m); SIG.append(s)

    MU = np.array(MU); SIG = np.array(SIG)
    if approx_x == False:
        return MU, SIG
    else:
        return MU, SIG, r_k

	


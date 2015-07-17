import numpy as np
from scipy.stats import norm
from scipy.misc import logsumexp

def kl_approx(mu1, sig1, mu2, sig2, size = 1):
	# compute KL divergence between 2 Gaussians
	if size == 1:
		KL = np.log(sig2) - np.log(sig1) - size + sig1 / sig2 + (mu1 - mu2) ** 2 / sig2
	else:
		is2 = np.inverse(sig2)
		KL = np.log(np.linalg.det(sig2)) - np.log(np.linalg.det(sig1)) - size + np.trace(np.dot(is2, sig1))
		KL += np.dot(np.dot((mu1 - mu2), is2), (mu1 - mu2))
	return KL / 2.0

def comp_log_likelihood(y, x, w, sig1, sig2, average = False):
	pdf1 = norm.pdf(y, loc = x, scale = np.sqrt(sig1))
	pdf2 = norm.pdf(y, loc = 0, scale = np.sqrt(sig2))
	Z = (1 - w) * pdf1 + w * pdf2
	if average == True:
		return np.log(Z).mean()
	else:
		return np.log(Z).sum()

def moment_numerical(mu, sig, y, sig1, sig2, w, average = False):
	length = 601
	a = 3 * np.sqrt(sig)
	inc = a / float(length - 1)
	x = np.linspace(-a + mu, mu + a, length)
	pdf = np.zeros(x.shape)
	for i in xrange(length):
		pdf[i] = comp_log_likelihood(y, x[i], w, sig1, sig2, average)
		pdf[i] += np.log(norm.pdf(x[i], loc = mu, scale = np.sqrt(sig)))
	logZ = logsumexp(pdf)
	pdf = np.exp(pdf - logZ)
	mean = (pdf * x).sum()
	cov = (x ** 2 * pdf).sum() - mean ** 2
	
	return cov, mean, logZ	
	#isx = 1 / cov - 1 / sig
	#misx = mean / cov - mu / sig
	#return isx, misx

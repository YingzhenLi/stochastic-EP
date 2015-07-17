import numpy as np
from comp_ep_functions import *

class mog(object):
	def __init__(self, size, prior_mu, prior_precision, std = 1.0, w = 0.5):
		self.size = size
		self.sig_noise = std ** 2
		if type(w) == float or (type(w) == list and len(w) == 1):
			self.w = np.array([w, 1 - w])
		else:
			self.w = np.array(w)
		self.w = self.w / float(self.w.sum())
		self.J = self.w.shape[0]	# number of groups
		if len(prior_precision.shape) > len(prior_mu.shape):
			self.full_cov = True
		else:
			self.full_cov = False
		self._init_params_prior(prior_mu, prior_precision)
		self._ep_param_initialsed = False

			
	def _init_params_prior(self, prior_mu, prior_precision):
		# params for prior
		self.pISx = prior_precision
		if self.full_cov is False:
			self.pMISx = prior_mu * self.pISx
		else:
			self.pMISx = prior_mu
			for j in xrange(self.J):
				self.pMISx[j] = np.dot(self.pISx[j], self.pMISx[j])
		
	def _init_ep_params(self, num_data, mode = 'full'):
		# local parameters for W
		self.num_data = num_data
		self.mode = mode
		# diagonal matrix first
		if self.full_cov is False:
			if mode == 'full':
				shape_is = [num_data, self.J, self.size]
				shape_mis = [num_data, self.J, self.size]
			else:
				shape_is = [self.J, self.size]
				shape_mis = [self.J, self.size]
		else:
			if mode == 'full':
				shape_is = [num_data, self.J, self.size, self.size]
				shape_mis = [num_data, self.J, self.size]
			else:
				shape_is = [self.J, self.size, self.size]
				shape_mis = [self.J, self.size]
		# local parameters
		self.isx = np.zeros(shape_is)
		self.misx = np.zeros(shape_mis)
		
		# global parameters
		if mode == 'full':
			self.ISx = self.isx.sum(0)
			self.MISx = self.misx.sum(0)
		else:
			self.ISx = self.isx * num_data
			self.MISx = self.misx * num_data
		self._ep_param_initialsed = True
		
	def check_positive_definiteness(self, matrix):
		if self.full_cov is False:
			success = (matrix > 0).all()
		else:
			success = not (np.any(np.linalg.eigvals(matrix) <= 0))
		return success
	
	def comp_current_moments(self):
		# return current moments
		if self.full_cov is False:
			Sn = 1.0 / (self.pISx + self.ISx)
			Mn = (self.pMISx + self.MISx) * Sn
		else:
			Sn = self.pISx + self.ISx
			Mn = (self.pMISx + self.MISx)
			for j in xrange(self.J):
				Sn[j] = np.linalg.inv(Sn[j])
				Mn[j] = np.dot(Sn[j], Mn[j])
		return Mn, Sn
		
	def comp_cavity(self, ind):
		if self.full_cov is False:
			Sn = 1.0 / (self.pISx + self.ISx - self.isx[ind])
			Mn = (self.pMISx + self.MISx - self.misx[ind]) * Sn
		else:
			Sn = self.pISx + self.ISx - self.isx[ind]
			Mn = (self.pMISx + self.MISx - self.misx[ind])
			for j in xrange(self.J):
				Sn[j] = np.linalg.inv(Sn[j])
				Mn[j] = np.dot(Sn[j], Mn[j])
		return Mn, Sn
	
	def comp_cavity_stochastic(self):
		if self.full_cov is False:
			Sn = 1.0 / (self.pISx + self.ISx - self.isx)
			Mn = (self.pMISx + self.MISx - self.misx) * Sn
		else:
			Sn = (self.pISx + self.ISx - self.isx)
			Mn = (self.pMISx + self.MISx - self.misx)
			for j in xrange(self.J):
				Sn[j] = np.linalg.inv(Sn[j])
				Mn[j] = np.dot(Sn[j], Mn[j])
		return Mn, Sn
		
	def update(self, ind, misx, isx, xi):
		tmp = self.ISx + xi * (isx - self.isx[ind])
		for j in xrange(self.J):
			if not self.check_positive_definiteness(tmp[j]):
				return 0
		self.ISx = self.ISx - self.isx[ind]
		self.MISx = self.MISx - self.misx[ind]
		self.isx[ind] = (1 - xi) * self.isx[ind] + xi * isx
		self.misx[ind] = (1 - xi) * self.misx[ind] + xi * misx
		self.ISx = self.ISx + self.isx[ind]
		self.MISx = self.MISx + self.misx[ind]
		
	def update_stochastic(self, minibatch, misx, isx, xi):
		tmp = self.ISx + xi * (isx - self.isx) * minibatch
		for j in xrange(self.J):
			if not self.check_positive_definiteness(tmp[j]):
				return 0
		self.ISx = self.ISx - self.isx * minibatch
		self.MISx = self.MISx - self.misx * minibatch
		self.isx = (1 - xi) * self.isx + xi * isx
		self.misx = (1 - xi) * self.misx + xi * misx
		self.ISx = self.ISx + self.isx * minibatch
		self.MISx = self.MISx + self.misx * minibatch
		self.isx = self.ISx / float(self.num_data)
		self.misx = self.MISx / float(self.num_data)
		
	def comp_local_updates(self, mu_t, sig_t, mu_c, sig_c):
		# compute local update
		# return no update if non-positive definite matrix
		
		# test1: diagonal matrix
		if self.full_cov is False:
			isx = 1.0 / sig_t - 1.0 / sig_c
			misx = mu_t / sig_t - mu_c / sig_c
			success = True
		else:
			isx = np.zeros(sig_c.shape)
			misx = np.zeros(mu_c.shape)
			for j in xrange(self.J):
				sig_t[j] = np.linalg.inv(sig_t[j])
				sig_c[j] = np.linalg.inv(sig_c[j])
				isx[j] = sig_t[j] - sig_c[j]
				misx[j] = np.dot(sig_t[j], mu_t[j]) - np.dot(sig_c[j], mu_c[j])
			success = True
				
		return misx, isx, success
		
	def comp_logZ(self, y):
		# TODO: need to be changed
		if self.mode == 'full':
			logZ = 0
			for ind in xrange(self.num_data):
				MU, SIG = self.comp_cavity(ind)
				pdf1 = norm.pdf(y[ind], loc = MU, scale = np.sqrt(SIG + self.sig1))
				pdf2 = norm.pdf(y[ind], loc = 0, scale = np.sqrt(self.sig2))
				logZ += np.log((1 - self.w) * pdf1 + self.w * pdf2)
		else:
			MU, SIG = self.comp_cavity_stochastic()
			pdf1 = norm.pdf(y, loc = MU, scale = np.sqrt(SIG + self.sig1))
			pdf2 = norm.pdf(y, loc = 0, scale = np.sqrt(self.sig2))
			logZ = np.log((1 - self.w) * pdf1 + self.w * pdf2).sum()
		return logZ
		
	def comp_gradient_logZ(self, y, delta = 0.001):
		SIG = 1 / (self.pISx + self.isx * (self.num_data - 1))
		MU = (self.pMISx + self.misx * (self.num_data - 1)) * SIG
		pdf1 = norm.pdf(y, loc = MU, scale = np.sqrt(SIG + self.sig1))
		pdf2 = norm.pdf(y, loc = 0, scale = np.sqrt(self.sig2))
		logZ1 = np.log((1 - self.w) * pdf1 + self.w * pdf2).sum()
		
		SIG = 1 / (self.pISx + (self.isx + delta) * (self.num_data - 1))
		MU = (self.pMISx + self.misx * (self.num_data - 1)) * SIG
		pdf1 = norm.pdf(y, loc = MU, scale = np.sqrt(SIG + self.sig1))
		pdf2 = norm.pdf(y, loc = 0, scale = np.sqrt(self.sig2))
		logZ2 = np.log((1 - self.w) * pdf1 + self.w * pdf2).sum()
		
		return (logZ2 - logZ1) / delta
		
	def train_ep(self, y, num_iter, learning_rate, mode):
		num_data = y.shape[0]
		# initialising ep parameters
		if self._ep_param_initialsed == False:
			self._init_ep_params(num_data, mode)
		# start training:
		if mode == 'adf':
			for epoch in xrange(num_iter):
				for ind in xrange(num_data):
					if self.full_cov is False:
						SIG = 1.0 / (self.pISx + self.ISx)
						MU = (self.pMISx + self.MISx) * SIG
					else:
						SIG = (self.pISx + self.ISx)
						MU = (self.pMISx + self.MISx)
						for j in xrange(self.J):
							SIG[j] = np.linalg.inv(SIG[j])
							MU[j] = np.dot(SIG[j], MU[j])
					mean, cov, r = gmm_updates([y[ind], self.w, self.sig_noise], SIG, MU, \
						approx_x = True, full_cov = self.full_cov)
					misx, isx, success = self.comp_local_updates(mean, cov, MU, SIG)
					if success:
						self.MISx += learning_rate * misx
						self.ISx += learning_rate * isx
		
		if mode == 'full':
			for epoch in xrange(num_iter):
				for ind in xrange(num_data):
					MU, SIG = self.comp_cavity(ind)
					mean, cov, r = gmm_updates([y[ind], self.w, self.sig_noise], SIG, MU, \
						approx_x = True, full_cov = self.full_cov)
					misx, isx, success = self.comp_local_updates(mean, cov, MU, SIG)
					if success:
						self.update(ind, misx, isx, learning_rate)
					
		if mode == 'stochastic':
			for epoch in xrange(num_iter):
				for i in xrange(num_data):
					ind = i#np.random.randint(num_data)
					MU, SIG = self.comp_cavity_stochastic()	
					mean, cov, r = gmm_updates([y[ind], self.w, self.sig_noise], SIG, MU, \
						approx_x = True, full_cov = self.full_cov)
					misx, isx, success = self.comp_local_updates(mean, cov, MU, SIG)
					if success:
						self.update_stochastic(1, misx, isx, learning_rate)
			
	def predict(self, X):
		# predict the cluster label
		SIG = (self.pISx + self.ISx)
		MU = self.pMISx + self.MISx
		if self.full_cov is False:
			SIG = 1.0 / SIG
			MU = MU * SIG
		else:
			for j in xrange(self.J):
				SIG[j] = np.linalg.inv(SIG[j])
				MU[j] = np.dot(SIG[j], MU[j])
		y_pred = np.zeros(X.shape[0], dtype = int)
		logZ_pred = np.zeros(X.shape[0])
		# TODO: need efficient implementation for processing multiple inputs together
		for i in xrange(X.shape[0]):
			r_k, logZ = gmm_updates([X[i], self.w, self.sig_noise], SIG, MU, pred = True, full_cov = self.full_cov)
			# take the max
			y_pred[i] = int(np.argmax(r_k))
			logZ_pred[i] = logZ
			
		return y_pred, logZ_pred
		

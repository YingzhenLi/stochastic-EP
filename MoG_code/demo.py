import numpy as np
from mog import mog
import sys, time
from comp_ground_truth import moment_numerical, kl_approx
from matplotlib import pyplot
from gmm import model as GMM
from try_bvn_ellipse import BVN, make_ellipses, gauss_ell
from sampling import gmm_sampling

color = ['r.', 'y.', 'g.', 'c.']

def kl_approx(mu1, sig1, mu2, sig2):
	# compute KL divergence between 2 Gaussians with full
	mu_diff = mu2 - mu1
	if len(sig2.shape) == 1:
	    KL = (sig1 / sig2).sum() - sig1.shape[0]
	    KL += (mu_diff ** 2 / sig2).sum()
	    KL += np.log(sig2).sum() - np.log(sig1).sum()
	else:
	    inv_sig2 = np.linalg.inv(sig2)
	    KL = np.trace(np.dot(inv_sig2, sig1)) - sig1.shape[0]
	    KL += np.dot(np.dot(mu_diff, inv_sig2), mu_diff)
	    KL += np.log(np.linalg.det(sig2)) - np.log(np.linalg.det(sig1)) 
	
	return KL / 2.0

def demo_clutter(seed, step, num_data, num_group, size, prior_precision, w, 
		std_noise, show=False, dim = [0, 1]):
	# generate data
	np.random.seed(seed*10)
	scale = 1.0
	MU_B = np.random.randn(num_group, size) * scale
	SIGMA_B = np.abs(np.random.randn(num_group, size))
	model = GMM(num_group, size, std_noise, MU_B, SIGMA_B, w)
    
    # Simulate_data
	print 'simulating training data...'
	data = model.simulate_data(num_data, seed=seed)
	X = data['X']	# observations
	y = data['y']	# cluster labels
	num_data = X.shape[0]
	
	# test data
	num_data_test = 1000
	print 'simulating test data...'
	data_test = model.simulate_data(num_data_test, seed=seed)
	X_test = data_test['X']
	y_test = data_test['y']
	
	if show:
		alpha = 0.55
		width = 5; hight = 3
		fig1, ax1 = pyplot.subplots(2, 2, figsize=(width, hight))
		pyplot.title('test')
		for n in xrange(X_test.shape[0]):
			ax1[0, 0].plot(X_test[n, dim[0]], X_test[n, dim[1]], \
				color[y_test[n]], alpha=alpha)
			ax1[0, 0].set_title('truth')
			ax1[0, 0].axis('off')
			
		# for training
		fig2, ax2 = pyplot.subplots(2, 2, figsize=(width, hight))
		pyplot.title('train')
		for n in xrange(X.shape[0]):
			ax2[0, 0].plot(X[n, dim[0]], X[n, dim[1]], color[y[n]], alpha=alpha)
			ax2[0, 0].set_title('truth')
			ax2[0, 0].axis('off')
		#pyplot.show()
	
	np.random.seed(0)
	prior_mu = np.random.randn(num_group, size)
	
	# sample from the true posterior
	level = 0.98
	sampling = True
	linewidth = 2.0
	if sampling:
		print 'computing the true posterior...'
		m_samp, var_samp, samp = gmm_sampling(X, y, num_group, prior_mu, \
			prior_precision, w, std_noise)
		# draw true posteior
		if show:
			bvn_full = BVN()
			cov_full = []; mean_full = []
			for k in xrange(num_group):
				cov_full.append(var_samp[k][np.ix_(dim, dim)])
				mean_full.append(m_samp[k][dim])			
			bvn_full.covars = cov_full
			bvn_full.mean = np.array(mean_full)
			make_ellipses(bvn_full, ax1[0, 0], level=level)
			for k in xrange(num_group):
				e1, e2 = gauss_ell(mean_full[k], cov_full[k], dim = [0, 1], \
					npoints = 200, level =level)
				ax1[0, 0].plot(e1, e2, 'k', linewidth=linewidth)
				ax2[0, 0].plot(e1, e2, 'k', linewidth=linewidth)
	
	# learning options
	mode = ['full', 'stochastic', 'adf']
	learning_rate = 0.1
	err = np.zeros([len(mode), num_track])	# error of posterior approximation
	ll = np.zeros([len(mode), num_track])	# for test likelihood
	t = np.zeros([len(mode), num_track])		# for training time
	
	# learning 'stochastic'
	name = {'full':'EP', 'stochastic':'SEP', 'adf':'ADF'}
	for m in xrange(len(mode)):
		print "fitting with %s..." % name[mode[m]]
		clutter_train = mog(size, prior_mu, prior_precision, std_noise, w)
		time_ep = time.time()
		for i in xrange(len(step)):
			clutter_train.train_ep(X, step[i], learning_rate, mode[m])
			t[m, i] = t[m, i] + time.time() - time_ep
			y_pred, logZ_pred = clutter_train.predict(X_test)
			y_pred_train, _ = clutter_train.predict(X)
			ll[m, i] = logZ_pred.mean()
			time_ep = time.time()
	
		if show:
			i = int(m >= 1); j = int(np.mod(m, 2) == 0)
			for n in xrange(X_test.shape[0]):
				ax1[i, j].plot(X_test[n, dim[0]], X_test[n, dim[1]], \
					color[y_pred[n]], alpha=alpha)
				ax1[i, j].set_title(name[mode[m]])
				ax1[i, j].axis('off')
				
			for n in xrange(X.shape[0]):
				current_color = y_pred_train[n]
				ax2[i, j].plot(X[n, dim[0]], X[n, dim[1]], \
					color[current_color], alpha=alpha)
				ax2[i, j].set_title(name[mode[m]])
				ax2[i, j].axis('off')
				
			bvn = BVN()
			cov = []; mean = []
			for k in xrange(num_group):
				if clutter_train.full_cov is True:
					cov.append(clutter_train.pISx[k][np.ix_(dim, dim)] \
						+ clutter_train.ISx[k][np.ix_(dim, dim)])
					cov[-1] = np.linalg.inv(cov[-1])
					mean.append(clutter_train.pMISx[k][dim] \
						+ clutter_train.MISx[k][dim])
					mean[-1] = np.dot(cov[-1], mean[-1])
				else:
					cov.append(clutter_train.pISx[k][dim] \
						+ clutter_train.ISx[k][dim])
					mean.append(clutter_train.pMISx[k][dim] \
						+ clutter_train.MISx[k][dim])
					mean[-1] = mean[-1] / cov[-1]
					cov[-1] = np.eye(2) / cov[-1]
					
			bvn.covars = cov
			bvn.mean = np.array(mean)
			make_ellipses(bvn, ax1[i, j], level=level)

			for k in xrange(num_group):
				e1, e2 = gauss_ell(mean[k], cov[k], dim = [0, 1], \
					npoints = 200, level =level)
				ax1[i, j].plot(e1, e2, 'k', linewidth=linewidth)
				ax2[i, j].plot(e1, e2, 'k', linewidth=linewidth)
	
	if show:
		pyplot.show()
		
	return ll, time_ep
		
if __name__ == '__main__':
	num_data = int(sys.argv[1])
	if len(sys.argv) > 2:
		size = int(sys.argv[2])
	else:
		size = 4
	if len(sys.argv) > 3:
		num_group = int(sys.argv[3])
	else:
		num_group = 4
	
	w = np.ones(num_group)
	if len(w) != num_group:
		w = np.ones(num_group)
	w = np.array(w); w = w / float(w.sum())
	std_noise = 0.5
	
	full_cov = True
	var_prior = 10.0
	if full_cov is False:
		prior_precision = np.ones([num_group, size]) / var_prior
	else:
		prior_precision = np.eye(size) / var_prior
		prior_precision = np.tile(prior_precision, (num_group, 1, 1))
	
	num_test = 1
	# for SEP paper figure, seed = [60]
	seed = [60]#np.arange(num_test) * 50
	num_mode = 3
	step = np.array([1, 2, 3, 4, 10, 30, 50, 100])
	num_track = len(step)
	ll = np.zeros([num_test, num_mode, num_track])
	time_ep = np.zeros([num_test, num_mode, num_track])
	dim = [0, 1]
	
	print 'settings:'
	print 'N_train_data = %d, dim = %d, N_clusters = %d, full Cov matrix = %s,' \
		% (num_data, size, num_group, full_cov)
	print 'total number of epochs = %d' % step.sum()
	
	for i in xrange(num_test):
		show = (i >= (num_test-1))
		ll[i], time_ep[i] = demo_clutter(seed[i], step, num_data, num_group, \
			size, prior_precision, w, std_noise, show=show, dim=dim)
	
	#print ll.mean(0)
	

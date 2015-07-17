import numpy as np
import pickle
from matplotlib import pyplot as plt

color = ['r', 'g', 'b']

def find_cluster(m1, m2, num_group):
	# find the corresponded cluster
	ind = np.array(range(num_group))
	res = []
	for i in xrange(num_group):
		mse = []
		for j in xrange(len(ind)):
			err = ((m1[i] - m2[ind[j]]) ** 2).sum()
			mse.append(err)
		mse = np.array(mse)
		index = np.argmin(mse)
		res.append(ind[index])
		ind = np.delete(ind, index)
	return res

def comp_error(size, num_group, num_test):
	filename = 'gmm_d%d_j%d_ntest%d' % (size, num_group, num_test)
	f = open(filename, 'r')
	option = ['EP', 'SEP', 'ADF']
	error_mean = {o:[] for o in option}
	error_var = {o:[] for o in option}
	time = 0
	
	# compute error
	for test in xrange(num_test):
		m_samp, var_samp, mean, cov, time_ep = pickle.load(f)
		for o in option:
			err_mean = []; err_var = []
			for i in xrange(mean[o].shape[0]):
				label = find_cluster(m_samp, mean[o][i], num_group)
				err_m = 0; err_v = 0
				for j in xrange(num_group):
					# taking averaged F-norm
					err_m += np.sqrt(((m_samp[j] - mean[o][i][label[j]]) ** 2).sum()) / float(num_group)
					err_v += np.sqrt(((var_samp[j] - cov[o][i][label[j]]) ** 2).sum()) / float(num_group)
				err_mean.append(err_m); err_var.append(err_v)
			err_mean = np.array(err_mean); err_var = np.array(err_var)
			error_mean[o].append(err_mean); error_var[o].append(err_var)
			time += time_ep
	
	# draw figure	
	fig, axs = plt.subplots(2, 1, figsize=(5, 4))
	fig.subplots_adjust(left=0.15, bottom=0.2)
	i = 0
	L = 1.5
	for o in option:
		error_mean[o] = np.array(error_mean[o])
		error_var[o] = np.array(error_var[o])
		
		mean_error_mean = error_mean[o].mean(0)
		var_error_mean = error_mean[o].var(0)
		mean_error_var = error_var[o].mean(0)
		var_error_var = error_var[o].var(0)
		
		xAxis = time[i] / float(num_test)
		axs[0].plot(xAxis, mean_error_mean, color = color[i], label=o, linewidth=L)
		axs[1].plot(xAxis, mean_error_var, color = color[i], label=o, linewidth=L)
		# plot variance
		plot_variance = False
		if plot_variance:
			axs[0].fill_between(xAxis, mean_error_mean + var_error_mean, \
				mean_error_mean - var_error_mean, color = color[i], alpha = 0.3)
			axs[1].fill_between(xAxis, mean_error_var + var_error_var, \
				mean_error_var - var_error_var, color = color[i], alpha = 0.3)
		i += 1
	
	#axs[0].set_title('error of mean')
	#axs[1].set_title('error of covariance')			
	axs[0].set_ylabel('averaged F-norm')
	axs[1].set_xlabel('time (s)')
	axs[1].set_ylabel('averaged F-norm')
	axs[0].set_yscale('log', basey=10)
	axs[1].set_yscale('log', basey=10)
	axs[0].set_xscale('log', basex=2)
	axs[1].set_xscale('log', basex=2)
	axs[0].get_xaxis().set_ticks([])
	axs[0].xaxis.set_ticks_position('bottom')
	axs[0].yaxis.set_ticks_position('left')
	axs[0].spines['top'].set_visible(False)
	axs[0].spines['right'].set_visible(False)
	axs[1].xaxis.set_ticks_position('bottom')
	axs[1].yaxis.set_ticks_position('left')
	axs[1].spines['top'].set_visible(False)
	axs[1].spines['right'].set_visible(False)
	lg = axs[0].legend(fontsize=12.0)
	lg.draw_frame(False)
	fig.frameon = False
	
	axs[0].set_ylim([0, 3])
    	
	plt.show()
    
if __name__ == '__main__':
	size = 4
	num_group = 4
	num_test = 5
	comp_error(size, num_group, num_test)
				

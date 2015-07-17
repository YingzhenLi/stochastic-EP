# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 21:49:03 2011

@author: josef
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


#function from David Cournapeau original scikits.learn em code
# To plot a confidence ellipse from multi-variate gaussian pdf
def gauss_ell(mu, va, dim = [0, 1], npoints = 100, level = 0.39):
    """ Given a mean and covariance for multi-variate
    gaussian, returns npoints points for the ellipse
    of confidence given by level (all points will be inside
    the ellipsoides with a probability equal to level)
    
    Returns the coordinate x and y of the ellipse"""
    
    c       = np.array(dim)

    if mu.size < 2:
        raise RuntimeError("this function only make sense for dimension 2 and more")

    if mu.size == va.size:
        mode    = 'diag'
    else:
        if va.ndim == 2:
            if va.shape[0] == va.shape[1]:
                mode    = 'full'
            else:
                raise DenError("variance not square")
        else:
            raise DenError("mean and variance are not dim conformant")

    # If X ~ N(mu, va), then [X` * va^(-1/2) * X] ~ Chi2
    chi22d  = stats.chi2(2)
    mahal   = np.sqrt(chi22d.ppf(level))
    
    # Generates a circle of npoints
    theta   = np.linspace(0, 2 * np.pi, npoints)
    circle  = mahal * np.array([np.cos(theta), np.sin(theta)])

    # Get the dimension which we are interested in:
    mu  = mu[dim]
    if mode == 'diag':
        va      = va[dim]
        elps    = np.outer(mu, np.ones(npoints))
        elps    += np.dot(np.diag(np.sqrt(va)), circle)
    elif mode == 'full':
        va  = va[c,:][:,c]
        # Method: compute the cholesky decomp of each cov matrix, that is
        # compute cova such as va = cova * cova' 
        # WARN: scipy is different than matlab here, as scipy computes a lower
        # triangular cholesky decomp: 
        #   - va = cova * cova' (scipy)
        #   - va = cova' * cova (matlab)
        # So take care when comparing results with matlab !
        cova    = np.linalg.cholesky(va)
        elps    = np.outer(mu, np.ones(npoints))
        elps    += np.dot(cova, circle)
    else:
        raise DenParam("var mode not recognized")

    return elps[0, :], elps[1, :]


#from sklearn gmm
def make_ellipses(bvn, ax, level=0.95):
    for n, color in enumerate('rgb'):
        v, w = np.linalg.eigh(bvn.covars[n][:2, :2])
        #print v, w
        u = w[0] / np.linalg.norm(w[0])
        angle = np.arctan(u[1]/u[0])
        angle = 180 * angle / np.pi # convert to degrees
        #v *= 39#25#9
        v = 2 * np.sqrt(v * stats.chi2.ppf(level, 2)) #JP
        ell = mpl.patches.Ellipse(bvn.mean[n, :2], v[0], v[1], 180 + angle,
                                  facecolor='none', 
                                  edgecolor=None, #color,
                                  ls='dashed',
                                  lw=3)
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(0.5)
        ax.add_artist(ell)
        
class BVN(object):
    pass
        

#cov1 = np.array([[1, 0.7],[0.7, 1]]) * 0.02
#mean1 = np.array([1,2]) /4.
#cov2 = np.array([[1, -0.3],[-0.3, 1]]) * 0.02
#mean2 = np.array([2,1]) /4.
#cov3 = np.array([[1, 0],[0, 0.5]]) * 0.02
#mean3 = np.array([2,2]) *2. /4.

#nobs = 1000
#rvs1 = np.random.multivariate_normal(mean1, cov1, size=nobs)
#rvs2 = np.random.multivariate_normal(mean2, cov2, size=nobs)
#rvs3 = np.random.multivariate_normal(mean3, cov3, size=nobs)
#print rvs1.mean(0)
#print rvs2.mean(0)
#print rvs3.mean(0)
#print np.cov(rvs1, rowvar=0)
#print np.cov(rvs2, rowvar=0)
#print np.cov(rvs3, rowvar=0)

#bvn = BVN()
#bvn.covars = [cov1, cov2, cov3]
#bvn.mean = np.vstack([mean1, mean2, mean3])



#fig = plt.figure()
#plt.plot(*rvs1.T, ls='.', color='r', alpha=0.25)
#plt.plot(*rvs2.T, ls='.', color='g', alpha=0.25)
#plt.plot(*rvs3.T, ls='.', color='b', alpha=0.25)
#level = 0.9
#for rvs, c in zip([rvs1, rvs2, rvs3], 'rgb'):
#    plt.plot(*rvs.T, ls='none', marker='.', color=c)#, alpha=0.25)
#ax = plt.gca()
#make_ellipses(bvn, ax, level=level)

#e1, e2 = gauss_ell(mean1, cov1, dim = [0, 1], npoints = 200, level =level)
#plt.plot(e1, e2, 'r')
#e1, e2 = gauss_ell(mean2, cov2, dim = [0, 1], npoints = 200, level = level)
#plt.plot(e1, e2, 'r')
#e1, e2 = gauss_ell(mean3, np.diag(cov3), dim = [0, 1], npoints = 200, 
#                   level = level)
#plt.plot(e1, e2, 'r')
#plt.show()

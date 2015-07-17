=========================================================
Stochastic EP for mixture of Gaussians (MoG) clustering
=========================================================
Author: Yingzhen Li

Comparing EP/ADF/SEP approximations of the posterior of the mean of 
Gaussian components. 
Details in the paper 
Yingzhen Li, Jose Miguel Hernandez-Lobato and Richard E. Turner.
"Stochastic Expectation Propagation". arXiv:1506.04132

Some of the codes are adjusted from:
ep-stan
https://github.com/gelman/ep-stan
try_bvn_ellipse.py
https://groups.google.com/d/msg/pystatsmodels/CmokUHssWiA/t1GOcI4dXkgJ

Required packages:
Stan
http://mc-stan.org/
pystan
http://pystan.readthedocs.org/en/latest/index.html
To run the demo, type
python demo.py num_train_data dim_data num_clusters

Example:
python demo.py 200 4 4


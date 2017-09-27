'''
Calculate the percentage of variance of expression when different 
number of genes are grouped together in each subgroup
written by Wenli Wu, last edited Sep 19, 2017
any questions related to the code, please email:wenliwu2018@u.northwestern.edu
'''

import numpy as np
import os,sys,glob
import pandas as pd
import macrogenomics

#### Set up data files ##############
cpmc = macrogenomics.cp_mc()
n_quantile_0 = np.exp(np.arange(1.5,7,0.2)).astype(int) # list of number of quantiles
r_s = [] # percentage of variance
C0,ld=cpmc.Ld_D()

#### Calculate the percentage of variance under each quantile ####
for n_quantile in n_quantile_0:
    Se = np.zeros((4,n_quantile))
    initial = np.zeros((4,n_quantile))
    for m in range(4):
        percentile_norm,Se[m,:] = cpmc.se_experiment(n_quantile,m)
    exp_g=np.exp(percentile_norm[1:])*cpmc.initial_aver
    y_fit=cpmc.se_g(exp_g) # model predicted result using g_function
    y_data = Se.mean(axis=0)[1:]*cpmc.D_fit/(C0*ld)
    r_s.append(cpmc.r_sq(y_data,y_fit))
r_s = np.array(r_s)
n_group = np.array(2445/n_quantile_0)
p = np.vstack((n_group,r_s)).transpose()
index = n_group<350
np.savetxt('percentOfVariance.csv',p[index,:],fmt='%3.2f,%.9f')

    
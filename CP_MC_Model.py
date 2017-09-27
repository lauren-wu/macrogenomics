'''
Calculate gene expression sensitivity based on the output of Monte Carlo 
and Brownian Dynamics output. The model predicted data is calculated 
using both direct simulation data fitting (sensitivity_model.csv) and 
the g function approximation(sensitivity_model_g.csv). The sentivity of expression 
measured by Microarray experiment is also calculated here as a validation
(sensitivity_experiment.csv).
written by Wenli Wu, last edited Sep 19, 2017
Any questions related to the code, please email:wenliwu2018@u.northwestern.edu
'''

import numpy as np
import os,sys,glob
import pandas as pd
import time
from sympy.solvers import solve
from sympy import Symbol
import sympy as sp
from scipy import stats
import scipy.special as sc
import macrogenomics
      
#######################################################################
#### Output the relation with initial expression ##############
#######################################################################
cpmc = macrogenomics.cp_mc()
n_quantile = 30# number of quantiles
Se_ex=np.zeros((4,n_quantile)) # sensitivity of microarray data
initial=np.zeros((4,n_quantile))

#### Output the experimental data ##########
for m in range(4):
    percentile_norm,Se_ex[m,:] = cpmc.se_experiment(n_quantile,m)
error = stats.sem(Se_ex)
C0,ld= cpmc.Ld_D()# C0=dD/dLd, calculate the prefactor to change sigma to D
p = np.vstack((percentile_norm[1:],Se_ex.mean(axis=0)[1:]*cpmc.D_fit/(C0*ld),error[1:])).transpose() # dln(E)/dD*D=(dx/dld*ld)*D/(C0*ld)
np.savetxt('sensitivity_experiment.csv',p,delimiter=',',fmt='%.9f')

#### Output the simulation prediction ##########
Se_model,expression,ratio=cpmc.se() # sensitivity predicted by the model
x = np.log(expression/cpmc.initial_aver)
y = Se_model
para2 = np.poly1d(np.polyfit(x,y,10))
x = np.linspace(-2.,2.8,100)
p = np.vstack((x,para2(x))).transpose()
np.savetxt("sensitivity_model.csv",p,delimiter=',',fmt='%.9f')

### Output function g ####
kappa,x,g_f =cpmc.g_fit(cpmc.D_fit) # calculate the critical point kappa
p = np.vstack((x,g_f)).transpose()
np.savetxt('g_function.csv',p,delimiter=',',fmt='%.9f')
print("the critical rate of expression of g function, kappa, is: ",kappa[0]*10**6," nM/s")

#### Output sensitivity calculated from g function
x = np.linspace(-2.,2.8,100)
exp_g = np.exp(x)*cpmc.initial_aver
y = cpmc.se_g(exp_g)
p = np.vstack((x,y)).transpose()
np.savetxt("sensitivity_model_g.csv",p,delimiter=',',fmt='%.9f')





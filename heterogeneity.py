'''
Calculate model predicted gene expression heterogeneity and the one measured from experiment
written by Wenli Wu, last edited Sep 19, 2017
any questions related to the code, please email:wenliwu2018@u.northwestern.edu
'''
import numpy as np
import os,sys,glob
import pandas as pd
from scipy import stats
import macrogenomics

#### Load the macrogenomic package #####################################
cpmc=macrogenomics.cp_mc()

########################################################################
#### Calculate the transcriptional heterogeneity predicted #############
#### by the model#######################################################
D_0=np.linspace(2.55,2.75,11) # fractal dimension
std_s=np.zeros(D_0.shape[0]) # intercellular trascriptional heterogeneity
i=0
for D in D_0:
    sigma=cpmc.phi_initial*(1-cpmc.phi_initial) # variance of continuous chromatin crowding distribution
    sigma_phi=sigma*(cpmc.Li/cpmc.r_min)**(D-3) # variance of the distribution of the local crowding condition in interaction volume
    Se,expression,ratio=cpmc.se()
    ASA=ratio**(-1/D) # 1e8 is the size of fractal
    std_s[i]=-ASA*cpmc.sec[2]*sigma_phi[2]
    i+=1

### Save the result
p=np.vstack((D_0,(std_s/std_s[0]))).transpose()
np.savetxt('heterogeneity_model.csv',p,delimiter=',',fmt='%.9f')

##### Calculate the Microarray data Heterogeneity #####################################
n_quantile=30 # number of quantiles
control_norm_mean=np.mean(cpmc.data[:,:4],axis=1) # get the microarray data for control
percentile_norm=np.percentile(control_norm_mean.flatten(),np.linspace(0,100,n_quantile)) # define the quantile by the expression of control (initial expression)
per_old=percentile_norm[4] # focus at the genes with initial expression around mean
n=0
for j in range(5,19):
    per=percentile_norm[j]   
    lis=((cpmc.data[:,:4]<=per) * (cpmc.data[:,:4]>=per_old))*1 # find the list of genes with expression around mean
    lis=lis.sum(axis=1)
    lis,=np.where(lis>1)
    std=np.zeros((lis.shape[0],4))
    # calculate the intercellular heterogeneity for each gene
    for s in range(lis.shape[0]):            
        h_function=np.zeros(cpmc.ld_exp.shape[0]) # the heterogeneity of each gene for 4 measurement groups with different D
        for i in range(cpmc.ld_exp.shape[0]):
            h_function[i]=np.var(cpmc.data[lis[s],i*4:i*4+4])
        a=np.vstack((cpmc.ld_exp,h_function)).transpose()
        a=a[a[:,0].argsort()]
        std[s]=a[:,1]/a[0,1]
    if n==0:
        H=std
    else:
        H=np.concatenate((H,std),axis=0)        
    n=n+1
    per_old=per
error=(stats.sem(H))/2/np.mean((H),axis=0)**(1/2)
D_measure=(a[:,0]-a[0,0])*1.6483+2.55 ## The measured D values in microarray experiment

## Save the result ###
p=np.vstack((D_measure,np.mean((H),axis=0)**(1/2),error)).transpose()
np.savetxt('heterogeneity_experiment.csv',p,delimiter=',',fmt='%.9f')


'''
Calculate Model predicted gene expression coefficient of variation
written by Wenli Wu, last edited Sep 19, 2017
any questions related to the code, please email:wenliwu2018@u.northwestern.edu
'''
import numpy as np
import os,sys,glob
import pandas as pd
from scipy import stats
import macrogenomics

## Simulation result
cpmc=macrogenomics.cp_mc()
sss=0
D_0=np.linspace(2.55,2.75,21)
for D in D_0:
	sss=sss+1
	pop,x,g_f=cpmc.g_fit(D)
	if sss==1:
		g_d=g_f	
	else:
		g_d=np.vstack((g_d,g_f))
n=2
cov_a=2**(1/2)*(1/g_d[:,n]-1)
p=np.vstack((D_0,cov_a/cov_a[0])).transpose()
np.savetxt('cov_model.csv',p,delimiter=',',fmt='%.9f')

##### Plot the Microarray data ##############
n_quantile=30## number of quantiles
initial=np.zeros(n_quantile)
control_norm_mean=np.mean(cpmc.data[:,:4],axis=1)
## define the quantile by the expression of control (initial expression)
percentile_norm=np.percentile(control_norm_mean.flatten(),np.linspace(0,100,n_quantile)) 
per_old=percentile_norm[4]
n=0

for j in range(5,18):
    per=percentile_norm[j]
    lis=((cpmc.data[:,:4]<=per) * (cpmc.data[:,:4]>=per_old))*1
    lis=lis.sum(axis=1)
    lis,=np.where(lis>1)
    cov=np.zeros((lis.shape[0],4))# coefficient of variation for genes with small range of initial expression
    for s in range(lis.shape[0]):
        var_each_gene=np.zeros(cpmc.ld_exp.shape[0])
        for i in range(cpmc.ld_exp.shape[0]):
            var_each_gene[i]=np.var(cpmc.data[lis[s],i*4:i*4+4])
        temp=np.vstack((cpmc.ld_exp,var_each_gene))
        a=temp.transpose()
        a=a[a[:,0].argsort()] # sort the var based on measured ld values            
        cov[s]=a[:,1]/a[0,1] 
    if n==0:
        text=cov
        mean_total=cpmc.data[lis,:]
    else:
        text=np.concatenate((text,cov),axis=0)       
        mean_total=np.concatenate((mean_total,cpmc.data[lis,:]),axis=0)
    n=n+1
    per_old=per
mean_temp=np.zeros(4)
for s in range(4):
    mean_temp[s]=np.mean(mean_total[:,s*4:s*4+4])
temp=np.vstack((cpmc.ld_exp,mean_temp))
a=temp.transpose()
a=a[a[:,0].argsort()]
mean_temp=a[:,1]
mean_temp=mean_temp/mean_temp[0]
error=(stats.sem(text))/2/np.mean((text),axis=0)**(1/2)/mean_temp
D=(a[:,0]-a[0,0])*1.6483+2.55
cov=np.mean((text),axis=0)**(1/2)/mean_temp

# save the output
p=np.vstack((D,cov,error)).transpose
np.savetxt('cov_experiment.csv',p,delimiter=',',fmt='%.9f')


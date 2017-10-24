''' 
Generate the mRNA expression ([mRNA]) as a function of phi under
different transcription factor concentration from the output of 
Brownian Dynamic simulation and Monte Carlo simulation based on 
the work done by Hiroaki et al. 2014

Output are the molecular feature (i.e. Transcription factor 
concentration, RNA polymerase concentration and binding site 
concentration. Saved as tot_con.csv as output), the maximum 
steady state mRNA concentration as a function of initial condition 
(max_mRNA_initial.csv), average crowding density phi at each molecular
condition (phi_initial.csv) and the second derivative of [mRNA] with 
respect to phi normalizedby the maximum mRNA concentration (second_derivative_TF_norm.csv).

written by Wenli Wu, last edited Sep 19, 2017
'''

import numpy as np
import matplotlib.pyplot as plt
import os,sys,glob
import pandas as pd
import time
from sympy.solvers import solve
from sympy import Symbol
import sympy as sp
import scipy.special as sc

## Enter the initial parameters ##
k_ns_t = 4.9e4 # mM-1s-1
k_ns_f = 3.6e4 # mM-1s-1
K_ns_D_TF_0 = 1 # mM
K_ns_D_RNAp_0 = 1 # mM
k_t = 0.0496e6 # mM-1s-1
k_f = 0.0324e6 # mM-1s-1
k_o = 0.991 # mM-1s-1
K_D_TF = 1e-6 # mM
K_D_RNAp = 1e-6 # mM
k_m = 0.001 # mM-1s-1
gamma = 8e-4 # s-1
mu = 3e-4 # s-1
TF_tot = 30e-6 # mM
RNAp_tot = 30e-6 #mM
O_tot = 30e-6 #mM
D_tot = 20 #mM


k_ns_o = K_ns_D_TF_0*k_ns_t
k_ns_b = K_ns_D_RNAp_0*k_ns_f
k_t_m = 500*(1e-6)**3*(8.7*1e-11/157.3*k_ns_o)**(1/2)
k_f_m = 500*(1e-6)**3*(8.7*1e-11/270.4*k_ns_b)**(1/2)
k_t = k_t_m*6.022*1e23
k_f = k_f_m*6.022*1e23
k_b = D_tot*K_D_RNAp/K_ns_D_RNAp_0*k_f
k_o = D_tot*K_D_TF/K_ns_D_TF_0*k_t

## Enter the simulation data ##
phi = np.arange(201)/201*0.505 # crowding density
f_TF_phi = 1 - 2.83*phi + 3.87*phi**2 - 4.11*phi**3;
f_RNAp_phi = 1 - 3.89*phi + 7.72*phi**2 - 7.72*phi**3;
f_cro_TF_phi = -3.2*phi - 2.0*phi**2;
f_cro_RNAp_phi = -3.7*phi - 2.7*phi**2;
f_cro_RNAps_phi = -2.6*phi - 4.6*phi**2;
f_ba_TF_phi = 2.5*phi**2;
f_ba_RNAp_phi = 3.1*phi**2;
f_ba_RNAps_phi = 0.1384*phi**2 + 9.203*phi**3;

## Enter simulation influenced parameters
K_ns_D_TF = K_ns_D_TF_0*np.exp(f_cro_TF_phi)
K_ns_D_RNAp = K_ns_D_RNAp_0*np.exp(f_cro_RNAp_phi)

k_ns_t_phi = k_ns_t*f_TF_phi*np.exp(-f_ba_TF_phi)
k_ns_f_phi = k_ns_f*f_RNAp_phi*np.exp(-f_ba_RNAp_phi)
k_ns_o_phi = k_ns_o*np.exp(f_cro_TF_phi)*f_TF_phi*np.exp(-f_ba_TF_phi)
k_ns_b_phi = k_ns_b*np.exp(f_cro_RNAp_phi)*f_RNAp_phi*np.exp(-f_ba_RNAp_phi)

k_t_phi = k_t*np.exp(f_cro_TF_phi/2)*f_TF_phi*np.exp(-f_ba_TF_phi/2)
k_f_phi = k_f*np.exp(f_cro_RNAp_phi/2)*f_RNAp_phi*np.exp(-f_ba_RNAp_phi/2)*np.exp(-f_ba_RNAps_phi)
k_o_phi = k_o*np.exp(f_cro_TF_phi/2)*f_TF_phi*np.exp(-f_ba_TF_phi/2)
k_b_phi = k_b*np.exp(f_cro_RNAps_phi)*np.exp(f_cro_RNAp_phi/2)*f_RNAp_phi*np.exp(-f_ba_RNAp_phi/2)*np.exp(-f_ba_RNAps_phi)

## Start the calculation of [mRNA](phi) ####
tot_con = 35*10**np.linspace(0,4,50)*1e-6 ## the intial molecular condition [mM]
mRNA_initial = np.zeros(len(tot_con)) # the [mRNA] at maximum under each iniyial condition at steady state[mM]
phi_range_real = np.zeros(len(tot_con))# the phi at mRNA_max under each initial condition [%]
sec = np.zeros(len(tot_con))# the second derivative of expression curve [mM]
mRNA_rate = np.zeros((len(tot_con),phi.size))# the mRNA expression under each initial condition [mM]

#### Start the transcription calcules ####################
n=0
C_I,C_II,O,TF_D,RNAp_D=sp.symbols('C_I,C_II,O,TF_D,RNAp_D')
for con in tot_con:
    TF_tot   = con
    RNAp_tot = con
    O_tot    = con
    mRNA     = np.zeros(phi.shape[0])
    ## Start the equation solving ##
    for i in range(phi.shape[0]):
        def F(C_I,C_II,O,TF_D,RNAp_D):
            'Find the steady state solution' 
            f1=O_tot-C_I-C_II-O
            f2=k_t_phi[i]*TF_D*O-k_o_phi[i]*C_I-k_f_phi[i]*RNAp_D*C_I+k_b_phi[i]*C_II
            f3=k_f_phi[i]*RNAp_D*C_I-(k_b_phi[i]+k_m)*C_II
            f4=(D_tot*k_ns_t_phi[i]/k_ns_o_phi[i]*(TF_tot-C_I-C_II)-k_m/k_ns_o_phi[i]*C_II)/(1+D_tot*k_ns_t_phi[i]/k_ns_o_phi[i])-TF_D
            f5=(D_tot*k_ns_f_phi[i]/k_ns_b_phi[i]*(RNAp_tot-C_II)-k_m/k_ns_b_phi[i]*C_II)/(1+D_tot*k_ns_f_phi[i]/k_ns_b_phi[i])-RNAp_D
            return(f1,f2,f3,f4,f5)
        result=sp.nsolve(F(C_I,C_II,O,TF_D,RNAp_D),(C_I,C_II,O,TF_D,RNAp_D),(1e-6,1e-6,1e-6,1e-6,1e-6))
        mRNA[i]=k_m/mu*float(result[1])
    #### Find the maximum phi ###################
    mRNA_rate[n,:]=mRNA
    index1=np.where(mRNA==max(mRNA))
    mRNA_initial[n]=mRNA[index1]       
    first_derivative = np.gradient(mRNA/mRNA[index1])/(phi[1]-phi[0])
    second_derivative = np.gradient(first_derivative)   
    sec[n]=second_derivative[index1]/(phi[1]-phi[0])
    phi_range_real[n]=phi[index1]
    #### End of change of phi ####################
    n=n+1
    
### Save the data ################################
np.savetxt('max_mRNA_initial_test.csv',mRNA_initial) # to generate mRNA production rate, the values in max_mRNA_initial.csv needs to be multiplied by mu, which is 3e-10 [s-1] in this case
np.savetxt('phi_initial_test.csv',phi_range_real)
np.savetxt('second_derivative_TF_norm_test.csv',sec)
np.savetxt('tot_con_test.csv',tot_con)


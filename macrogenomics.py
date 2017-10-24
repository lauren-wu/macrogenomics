'''
Calculate the gene expression sensitivity the scaling of chromatin density model
written by Wenli Wu, last edited Sep 19, 2017
any questions related to the code, please send the email to : wenliwu2018@u.northwestern.edu

'''
import numpy as np
import os,sys,glob
import pandas as pd
from sympy import Symbol
import sympy as sp
from scipy import stats
import scipy.special as sc
from scipy.optimize import curve_fit

class cp_mc:
    '''
    The Chromatin Packing Molecular Crowding Model
	''' 
    ## load experiment Data
    data = pd.read_excel('full_genes.xlsx',header=None).values
    data_norm = np.log(data/np.mean(data))
    ld_exp = np.array([1.001, 1., 0.9933, 0.9151]) # Relative Ld value
 
    ## Molecular Crowding output and Chromatin Packing parameters
    sec          =np.genfromtxt('second_derivative_TF_norm_test.csv')
    mRNA_initial =np.genfromtxt('max_mRNA_initial_test.csv')
    phi_initial  =np.genfromtxt('phi_initial_test.csv')
    tot_con      =np.genfromtxt('tot_con_test.csv')
    D_fit        =2.68 ## The average D of cells
    Li0          =15e-9 # Size of interaction volume for single base pair
    r_max        =1000e-9 # maximum radius of chromatin packing scale
    r_min        =1e-9 # size of elementary particle
    len_gene     =6e3 # average length of genes
    initial_aver =0.00056 # the average initial expression rate
    Li           =(Li0+len_gene**(1/D_fit)*r_min) # the total interaction volume
    
    
    def Ld_D(self):
        '''
        the PWS code to transfer sigma(ld) to D
        output:
        C0: the derivative of D as respect to Sigma (ld) 
        Ld: PWS measurement parameter Ld(sigma)
        '''
        ################# Lusik's code ########################
        Ln  = 1.50e-6# microns
        D   = 2.68 # D of the correlation function
        L   = 3.0e-6  # sample thickness in microns
        ##### written by L. Cherkezyan, last edited Jan 23 2016
        ##### Inplemented in Python by Wenli Wu
        #######################################################

        pi      = 3.14
        sigma_n = 0.05/1.53
        # parameters of PWS system:
        n1      = 1
        n2      = 1.53
        NA      = 0.6
        G_top   = abs(((n1-n2)/(n1+n2))**2)# Fresenel reflection at top surface
        G       = abs((n1-n2)*2*n1*2*n2/((n1+n2)**3))/G_top
        k1      = n2*2*pi/0.7e-6# max wavelength = 0.7um
        k2      = n2*2*pi/0.5e-6# max wavelength = 0.5um
        kc      = (k1+k2)/2 # central wavenumber
        r_min = 1*1e-9 #um, define the minimum length scale at which RI exists
        dlist = np.zeros(100) 
        ldList = np.zeros(100)
        D_list=np.linspace(2,2.99,100)
        qqq=0
        for D in D_list:
            An = sigma_n**2*(r_min/Ln)**(1.5-D/2)/sc.kv(D/2-1.5,r_min/Ln)# this definition of An keeps RI variance constant independent of Ln and D
            x  = kc*Ln
            Sigma_R_2 = 2*An*(G**2*kc*L)/np.sqrt(pi)*sc.gamma(D/2)*x/((D-2)*2**(1.5-D/2))*((1+4*x**2)**(1-D/2)-(1+x**2*(4+NA**2))**(1-D/2))
            Sigma_L_2 = An*2**((D-9)/2)*G**2*sc.gamma((D-3)/2)*(1-(1+x**2*NA**2)**(1.5-D/2))
            Sigma  = np.sqrt(Sigma_R_2+Sigma_L_2)
            Ld = Sigma/np.sqrt(kc*L/n2)## sigma is variance of spectrum or mass??
            ldList[qqq] = Ld
            qqq = qqq+1
        num = np.argmin(abs(D_list-self.D_fit))
        L_d = ldList[num] ## Ld for state one is 1 micron
        dL_d = ldList[num+1]-ldList[num]
        dD = D_list[num+1]-D_list[num]
        C0 = dD/dL_d
        return(C0,L_d)
    def se(self):
        '''
        Calculate the gene expression sensitivity using cp-mc model
        output:
        Se: sensitivity as a function of initial average expression
        expression: initial average expression in each subgroups (decided by the molecular condition)
        ratio: the ratio between the total mass of fractal and the mass of elementary particle, M_f/M_in
        '''
        ratio = (self.r_max/self.r_min)**self.D_fit # M_f/M_min
        sigma = self.phi_initial*(1-self.phi_initial) # variance of crowding density at each point
        sigma_phi = sigma*(self.Li/self.r_min)**(self.D_fit-3) # variance of average crowding density in each interaction volume 
        d_sigma_phi = sigma_phi*(np.log(self.Li/self.r_min)+(3-self.D_fit)/self.D_fit**2*self.len_gene**(1/self.D_fit)*np.log(self.len_gene)*(self.r_min/self.Li)) # the derivative of sigma_phi as respect to D
        Se = (1/self.D_fit**2*np.log(ratio)+1/2*self.sec/(1+1/2*self.sec*sigma_phi)*d_sigma_phi)*self.D_fit # gene transcription Sensitivity
        initial_expression = self.initial_aver 
        expression = self.mRNA_initial*(1+1/2*self.sec*sigma_phi) # the average expression of gene share similar molecular properties
        return (Se,expression,ratio)

    def se_experiment(self,n_quantile,m):
        '''
        calculate sensitivity data from experiment
        input:
        n_quantile: number of quantile we want to divide the data into
        m: the number of treatment group. In this microarray study, m=1 is control, m=2,3,4 are treated groups
        output:
        percentile_norm: the expression at each quantile
        Se: sensitivity of expression with each quantile
        '''
        ## define the quantile by the expression of control (initial expression)
        control_norm_mean = self.data_norm[:,m] #the log normalized expression for control group
        percentile_norm = np.percentile(control_norm_mean.flatten(),np.linspace(0,100,n_quantile+1)[1:]) # devide the expression to subgroups based on the quantile of control
        per_old = -100
        Se = np.zeros(n_quantile) # sensitivity of expression for experiment data
        
        # find the sensitivity of expression in each quantile
        for j in range(n_quantile):
            per = percentile_norm[j]
            lis, = np.where((control_norm_mean <= per)&(control_norm_mean>=per_old))         
            Se_temp = np.zeros(lis.shape[0])
            for s in range(lis.shape[0]):
                ld_function = np.zeros(self.ld_exp.shape[0])
                for i in range(self.ld_exp.shape[0]):
                    ld_function[i] = self.data_norm[lis[s],(i)*4+m]
                temp = np.vstack((self.ld_exp,ld_function))
                a = temp.transpose()
                a = a[a[:,0].argsort()]
                Se_temp[s] = np.polyfit(a[:,0],a[:,1],1)[0]
            Se[j] = np.mean(Se_temp)
            per_old = per
        return (percentile_norm,Se)

    def func2(self,x,a):
        '''
        the fitting function for g
        input:
        x: intial mRNA expression at average phi (~ the maximum phi), epsilon_0
        a: the fitting parameter of g function, kappa=a^2*3e-4, where 3e-4 is the degradation rate of mRNA,
        used to change steady state expression concentration to transcription rate 
        output: the average expression of genes sharing similar molecular properties, epsilon_bar
        ''' 
        return a/x**(1/2)

    def g_function(self,x,popt,sigma):
        '''
        simplied format of g function
        input:
        x: the average expression of genes sharing similar molecular properties, epsilon_bar
        popt: a, where kappa=a^2*3e-4
        sigma: the variance of average crowding density within each interaction volume
        output: the g function as a function of x (epsilon_bar)
        '''
        return 8*x/(8*x+popt**2*sigma**2+(16*popt**2*x*sigma**2+popt**4*sigma**4)**(1/2))
       
    def g_fit(self,D):
        '''
        g function relating epsilon_bar and epsilon(m,phi_bar)
        input:
        D: fractal dimension
        output:
        popt: fitting parameter kappa, the critical expression rate
        x: the average expression of genes sharing similar molecular properties, epsilon_bar
        g_f: the fitted g function as respect to x
        '''
        epsilon_0=self.mRNA_initial # the range of mRNA_initial that makes epsilon_bar positive
        ydata=self.sec
        popt, pcov = curve_fit(self.func2, epsilon_0, ydata) # popt is the square root of kappa devided by the degradation rate of mRNA
        sigma=self.phi_initial*(1-self.phi_initial)*(self.Li/self.r_min)**(D-3)
        epsilon_bar=self.mRNA_initial*(1+1/2*sigma*self.sec)
        x=epsilon_bar
        g_f=self.g_function(x,popt[0],sigma)
        popt=popt**2*3e-4
        return (popt,x,g_f) # kappa=popt^2*degradation rate of mRNA
        
    def se_g(self,x):
        '''
        The sensitivity curve calculated by the approximation of g function
        input:
        x: the average expression of genes sharing similar molecular properties, epsilon_bar
        output:
        se: the sensitivity calculated by the g function
        '''
        ratio = (self.r_max/self.r_min)**self.D_fit
        pop,a,b=self.g_fit(self.D_fit) # fit the g function to find critical point kappa
        pop=(pop/3e-4)**(1/2)
        sigma=self.phi_initial[25]*(1-self.phi_initial[25])*(self.Li/self.r_min)**(self.D_fit-3)      
        se=(1-1/self.g_function(x,pop[0],sigma))*(self.D_fit*np.log(self.Li/self.r_min)+(3-self.D_fit)/self.D_fit*(self.r_min/self.Li)*self.len_gene**(1/self.D_fit)*np.log(self.len_gene))+1/self.D_fit*np.log(ratio)
        return se

    def r_sq(self,y_data,y_fit):
        '''
        the Percentage of variance for nonlinear fit in the model
        '''
        ssr = np.sum((y_fit-y_data-np.mean(y_fit-y_data))**2)
        rrio = np.sum((y_data-np.mean(y_data))**2)
        return 1-ssr/rrio





# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 13:45:02 2020

@author: jlamon02
"""

import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr

def MSE(S,O):
    MSE_hat = np.mean((S-O)**2)
    return MSE_hat

def NSE(S,O):
    NSE_hat = 1.-MSE(S,O)/np.std(O,ddof=1)
    return NSE_hat

def LNSE(S,O):
    V = np.log(S)
    U = np.log(O)
    LNSE_hat = NSE(U,V)
    return LNSE_hat

def KGE_prime(S,O):
    beta = np.mean(S)/np.mean(O)
    alpha = np.std(S,ddof=1)/np.std(O,ddof=1)
    rho = pearsonr(S,O)[0]
    KGE_hat = 1.-np.sqrt((beta-1.)**2+(alpha-1)**2+(rho-1)**2)
    return KGE_hat

def PVSE_prime(S,O):
    S = np.sort(S)
    O = np.sort(O)
    beta = np.mean(S)/np.mean(O)
    alpha = 1 - 0.5*np.mean(np.abs(S/np.mean(S)-O/np.mean(O)))
    rho,p = spearmanr(S,O)
    PVSE_hat = 1.-np.sqrt((beta-1.)**2.+(alpha-1.)**2.+(rho-1.)**2.)
    return PVSE_hat

def LBE(S,O,prime=False):
    mu_O,var_O,mu_S,var_S,rho = BLN3_moments(O,S)
    
    Co = np.sqrt(var_O)/mu_O#Eq. 18
    alpha = np.sqrt(var_S)/np.sqrt(var_O)
    delta = 1.-mu_S/mu_O
    beta = 1.-delta

    if prime == False:
        LBE_hat = 2.*alpha*rho-alpha**2.-(delta**2)/(Co**2) #Lamontagne et al. (2020) eq. 18a
    elif prime == True:
        LBE_hat = 1.-np.sqrt((beta-1.)**2.+(alpha-1.)**2.+(rho-1.)**2.) #Lamontagne et al. (2020) eq. 18b
    return LBE_hat

def LBEm(S,O,period,prime=False):
    mu_Om,mu_Sm,var_Om,var_Sm,mu_SOm = mixture_moments(S,O,period,prime=False)
    alpha_m = np.sqrt(var_Sm)/np.sqrt(var_Om)#Lamontagne et al. 2020 eqn. 20
    delta_m = 1.-mu_Sm/mu_Om
    Co_m = np.sqrt(var_Om)/mu_Om
    r_m = (mu_SOm-mu_Om*mu_Sm)/(np.sqrt(var_Om)*np.sqrt(var_Sm))
    
    if prime == False:
        LBEm_hat = 2*alpha_m*r_m-alpha_m**2-(delta_m**2)/(Co_m**2)#Lamontagne et al. 2020 eqn. 20a
    elif prime == True:
        LBEm_hat = 1-np.sqrt((delta_m**2)+(alpha_m-1)**2+(r_m-1)**2)#Lamontagne et al. 2020 eqn. 20b
    return LBEm_hat

#Supporting code
def BLN3_moments(O,S):
    U,tau_O = LN3_trans(O)
    V,tau_S = LN3_trans(S)
    
    U_bar = np.mean(U)
    V_bar = np.mean(V)
    var_U = np.var(U,ddof=1)
    var_V = np.var(V,ddof=1)
    cov_UV = np.cov(U,V,rowvar=False,ddof=None)[0,1]
    
    mu_O = tau_O+np.exp(U_bar+var_U/2.)
    mu_S = tau_S+np.exp(V_bar+var_V/2.)
    var_O = np.exp(2.*U_bar+var_U)*(np.exp(var_U)-1.)
    var_S = np.exp(2.*V_bar+var_V)*(np.exp(var_V)-1.)
    
    rho = (np.exp(cov_UV)-1.)/np.sqrt((np.exp(var_U)-1.)*(np.exp(var_V)-1.))
    return mu_O,var_O,mu_S,var_S,rho

def stedinger_tau(X):
    if X.min()+X.max()-2.*np.median(X)>0:
        tau = (X.min()*X.max()-np.median(X)**2)/(X.min()+X.max()-2.*np.median(X))
    else:
        tau = 0.
    return tau

def LN3_trans(X):
    tau = stedinger_tau(X)
    Y = np.log(X-tau)
    return Y,tau

def mixture_moments(S,O,period,prime=False):
    m = np.size(np.unique(period))
    mu_O = np.zeros([m,1])
    var_O = np.zeros([m,1])
    mu_S = np.zeros([m,1])
    var_S = np.zeros([m,1])
    rho = np.zeros([m,1])
    
    for i in range(m):
        mu_O[i],var_O[i],mu_S[i],var_S[i],rho[i] = BLN3_moments(O[period==i],S[period==i])
    mu_O_hat = np.mean(mu_O)
    mu_S_hat = np.mean(mu_S)
    
    var_O_hat = np.mean(var_O+mu_O**2)-mu_O_hat**2
    var_S_hat = np.mean(var_S+mu_S**2)-mu_S_hat**2
    
    mu_SO_hat = np.mean(mu_S*mu_O+rho*np.sqrt(var_S)*np.sqrt(var_O))
    return mu_O_hat,mu_S_hat,var_O_hat,var_S_hat,mu_SO_hat


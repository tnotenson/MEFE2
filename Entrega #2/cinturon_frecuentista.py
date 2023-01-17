#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 20:45:02 2021

Herramienta para la construcción de cinturones de confianza frecuentista

@author: tomasnotenson
"""

# Titulo: Herramienta para la construcción de cinturones de confianza frecuentista

# Parte A
#  1. Generar n números aleatorios xi con distribución exponencial con un
#  valor dado de τ y calcular un estadístico t igual al promedio de los xi.
#  2. Repetir el punto anterior, un numero NExperimentos de veces y guardar
#  el valor del estadístico t obtenido cada vez en un histograma.
#  3. Para ese histograma, encontrar un valor de tmin y tmax tal que entre
#  ellos se encuentre una fracción CL del total de los eventos.
#  4. Repetir los tres puntos anteriores para 100 valores de tau entre 0 y 10
#  y graficar  el cinturón de confianza en el plano τ vs t.

# Parte B
#  Utilizar la herramienta generada en la Parte A para otros dos estadísticos:
#  1. la mediana 
#  2. t=Σi  (xi+xi**2+xi**3+xi**4)/n

# Parte C
# Discutir que pasa si como "estadístico" t utilizamos el q de Wilks visto en clase. 

# Parte D
# Utilizar Wilks para calcular el intervalo de confianza para el parámetro τ de la  exponencial con un dado CL.

# Parte E
# (Opcional) Calcular la cobertura de los intervalos calculados en los ítems anteriores en función del parámetro τ 

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from statistics import median
from time import time

plt.rcParams.update({
"text.usetex": True,
"font.family": "sans-serif",
"font.sans-serif": ["Helvetica"]})

def freedman_diaconis(data, returnas="width"):
    """
    Use Freedman Diaconis rule to compute optimal histogram bin width. 
    ``returnas`` can be one of "width" or "bins", indicating whether
    the bin width or number of bins should be returned respectively. 


    Parameters
    ----------
    data: np.ndarray
        One-dimensional array.

    returnas: {"width", "bins"}
        If "width", return the estimated width for each histogram bin. 
        If "bins", return the number of bins suggested by rule.
    """
    data = np.asarray(data, dtype=np.float_)
    IQR  = stats.iqr(data, rng=(25, 75), scale="raw", nan_policy="omit")
    N    = data.size
    bw   = (2 * IQR) / np.power(N, 1/3)

    if returnas=="width":
        result = bw
    else:
        datmin, datmax = [data.min(), data.max()]
        datrng = datmax - datmin
        # print(datrng)
        result = int((datrng / bw) + 1)
    return(result)

def prom(x):
    return np.mean(x)

def falopa(x):
    return np.sum([xi**1+xi**2+xi**3+xi**4 for xi in x])/len(x)

def LL(tau,x):
    # print(tau)
    return np.prod(1/tau*np.exp(-x/tau))

def qWilks(tau,x):
    tauMLE = prom(x)
    q = -2*np.log(LL(tau,x)/LL(tauMLE,x))
    return q

def cinturon_frecuentista_exp(par,funt,n,Nexp,CL):
        
    CI_t = np.zeros((len(par), 2))
    
    start_time = time()
    # print(type(par))
    for i,tau in enumerate(par):
            
        # A)2)
        ts = np.zeros((Nexp))
        for iexp in range(Nexp):
            # A)1)
            xs = np.random.exponential(tau, n)
            if funt == qWilks:
                t = funt(tau,xs)
            else:
                t = funt(xs)
            # print(t)
            ts[iexp] = t
        # print(ts)
        # nbines = freedman_diaconis(ts, 'bins') # calculo la cantidad de bines óptima segun F-D
        # ft, bines = np.histogram(ts, density=True, bins=nbines)[:2] # defino los parámetros para el histograma
        # w = np.diff(bines) # Ancho de los bins
        # # #######################################################################
        # # grafico para chequear cada paso
        # plt.figure()
        # plt.bar(bines[:-1], ft, width=w, align='edge', ec='k') # histograma
        # plt.xlabel(r'$t$')
        # plt.grid(alpha=0.5)
        # plt.savefig('histograma_t.png', dpi=80)
        #######################################################################
        
        # ft_normed = ft/(w*np.sum(ft))
        
        # A)3)
        # CL = 0.683
        CL_l = (1-CL)/2
        CL_u = 1-CL_l
        cuantiles = np.quantile(ts, q=[CL_l, CL_u])
        CI_t[i,:] = cuantiles
        #######################################################################
        # grafico para chequear
        # plt.vlines(cuantiles[0], min(ft), max(ft), color='r'); plt.vlines(cuantiles[1], min(ft), max(ft), color='r')
        #######################################################################
        
        #######################################################################
        # chequeo el CL
        # arr1 = np.where(bines <= cuantiles[1], bines, 0*bines)
        # arr2 = np.where(cuantiles[0] <= bines, bines, 0*bines)
        # indx_integral = np.where(arr1==arr2)
        # bines_integral = bines[indx_integral]
        # ft_integral = ft_normed[indx_integral]
        
        # CL_real = np.sum(w[0]*ft_integral)
        # print(f'Real CL = {CL_real}')
        
        # t_l, t_u = [min(bines_integral), max(bines_integral)]
        #######################################################################
        print(f't confidence interval: [{cuantiles[0]}; {cuantiles[1]}]')
    final_time = time()
    print(f'Total duration "cinturon frecuentista exp" A): {final_time - start_time} seg')
    return CI_t
#%%
# defino parámetros, cantidad de datos y de experimentos
taus = np.linspace(0.001,10,100)
n = 100
Nexp = 1000
CL = 0.683

# elijo el estadístico a utilizar
# A)
funt = prom 
# B)
# funt = median
# funt = falopa 
# C)
# funt = qWilks
if funt == prom:
    sfun = 'prom'
elif funt == median:
    sfun = 'median'
elif funt == falopa:
    sfun = 'falopa'
elif funt == qWilks:
    sfun = 'qWilks'
    
CI_t = cinturon_frecuentista_exp(taus,funt,n,Nexp,CL)

#%% Chequeo y ploteo
tau = 10
xs = np.random.exponential(tau, n)
tp = funt(xs)
##############################################################################
# chequeo el CL en el límite asintótico
# tvar = tp**2/n
# tsd = np.sqrt(tvar)
# t_dist= stats.norm()#gamma(n,n*tau)

# CL_l = (1-CL)/2
# CL_u = 1-CL_l

# a, b = [t_dist.ppf(CL_l), t_dist.ppf(CL_u)]
# CI = [tp + a*tsd, tp+b*tsd]
# print(f'Analytic confidence interval: [{CI[0]};{CI[1]}]')
# print(f'AsyLim confidence interval: [{CI_t[20][0]};{CI_t[20][1]}]')
##############################################################################
#%% ploteo t vs τ y τ vs t en la misma figura
fig, ax = plt.subplots(2)
ax[0].plot(taus, CI_t[:,0],'.b')
ax[0].plot(taus, CI_t[:,1],'.b')
# ax[0].hlines(tp, min(taus), max(taus), 'r')
ax[0].set_xlabel(r'$\tau$')
ax[0].set_ylabel(r'$t$')
ax[0].grid()

ax[1].plot(CI_t[:,0], taus,'.b')
ax[1].plot(CI_t[:,1], taus, '.b')
# ax[1].vlines(tp, min(taus), max(taus), 'r')
ax[1].set_xlabel(r'$t$')
ax[1].set_ylabel(r'$\tau$')
ax[1].grid()
plt.savefig(f't_vs_tau_n{n}_CL{CL}_Nexp{Nexp}_fun{sfun}.png', dpi=80)
#%% ploteo solamente τ vs t
plt.plot(CI_t[:,0], taus,'.b')
plt.plot(CI_t[:,1], taus, '.b')
# plt.vlines(tp, min(taus), max(taus), 'r')
plt.xlabel(r'$t$')
plt.ylabel(r'$\tau$')
plt.grid()
plt.savefig(f't_vs_tau_n{n}_CL{CL}_Nexp{Nexp}_fun{sfun}.png', dpi=80)
plt.close()
#%%

# D) Intervalo usando la verosimilitud: LL(ci_limits) = LLmax-1/2 ///////////

def LLexp(x,tau):
    tauMLE = prom(x)
    n = len(x)
    return -n*tauMLE/tau - n*np.log(tau)

n = 100 
tau = 1
xs = np.random.exponential(tau,n)
tauMLE = prom(xs)
dom = np.linspace(0.1, 5*tau,10000)
img = LLexp(xs,dom); y = LLexp(xs,tauMLE) - 1/2
arr1 = np.where(img >= y)
print(arr1)
xmin = dom[arr1[0][0]]; xmax = dom[arr1[0][-1]]

plt.figure()
plt.plot(dom, img, '-b', tauMLE, LLexp(xs,tauMLE), '.r')
plt.hlines(y,min(dom),max(dom), ls='dashed', lw = 0.8, alpha=0.8)
plt.vlines(xmin,min(img),max(img), ls='dashed', lw = 0.8, alpha=0.8)
plt.vlines(xmax,min(img),max(img), ls='dashed', lw = 0.8, alpha=0.8)
plt.fill_between(dom, min(img), max(img), where=img >= y,
                color='blue', alpha=0.2)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$log(\mathcal{L})$')
plt.grid()
plt.show()
plt.savefig(f'LL_vs_tau_n{n}_Wilks.png', dpi=80)
#%%
# def IsInside2(x,mu, CI=False):
#     # Esta función calcula el intervalo de confianza de Wilks
#     # para obtener el intervalo de confianza de Wilks, se debe cambiar el
#     # parámetro CI a True y tomar el segundo parámetro
#     dom = np.linspace(0.001, 10*mu,10000)
#     img = LLvec(x,dom); y = max(img) - 1/2
#     arr1 = np.where(img >= y)
#     xmin = dom[arr1[0][0]]; xmax = dom[arr1[0][-1]]
#     if not CI:
#         return (xmin <= mu and mu <= xmax)
#     else:
#         return [(xmin <= mu and mu <= xmax), [xmin,xmax]]
#%%
# Este script calcula el intervalo frecuentista para el mu de una variable aleatoria de Poisson.

def cobertura_vs_par(par, IsInside=IsInside2, funt=prom):
    
    start_time = time()
    
    probs = np.zeros((nscan_points))
    # Bucle sobre mu, desde mu_min a mu_max.
    # Para cada mu calcula la cobertura del intervalo y lo guarda en probs
    for i in range(0,nscan_points):
        mu = par[i]
        # print(mu)
       
        Nmax = int(mu+10*np.sqrt(mu))+1 # esperanza mas 10 veces sigma
        xs = np.random.exponential(mu, n)
        # nbines = freedman_diaconis(xs, 'bins') # calculo la cantidad de bines óptima segun F-D
        # ft, bines = np.histogram(xs, density=True, bins=nbines)[:2] # defino los parámetros para el histograma
        # w = np.diff(bines) # Ancho de los bins
        
        probInside2 = 0
        for j in range(nbines):
            x = bines[j]
            # Barro desde n observado hasta la esperanza mas 10 veces sigma    
            if IsInside == IsInside2:
                inside2 = IsInside(x,mu)
            elif IsInside == IsInside_CI:
                inside2 = IsInside(x,mu,funt)
                
            prob = stats.expon.pdf(x,mu);

    # Si nobs esta dentro del intervalo, agrega la probabilidd de ese caso particular, i.e. Poisson(nobs,mu).
        if (inside2): probInside2 += prob*w;
        
        probs[i] = probInside2
    final_time = time()
    print(f'TOTAL duration D: {final_time - start_time} seg')
    return probs
#%% Tomo un ejemplo del CI de Wilks
Be, CI = IsInside2(tp, tau, CI=True)
print(Be, CI)
#%% Cobertura vs τ
mu_min = 0.01
mu_max = 6.0
nscan_points= 1000
mus = np.linspace(mu_min, mu_max,nscan_points)
CL = 0.6827

probs = cobertura_vs_par(mus, funt=funt)    

plt.figure()    
plt.plot(mus,probs,'-r')
plt.hlines(CL, min(mus), max(mus), ls='dashed', lw = 0.8, alpha=0.8)
plt.xlabel(r'$\mu$')
plt.ylabel(r'$CL$')
plt.grid()
plt.savefig(f'cobertura_vs_tau_mumin{mu_min}_mumax{mu_max}_scanpoints{nscan_points}_fun{sfun}.png', dpi=200)
#%% E) ####IN PROCESS##### 
from scipy.interpolate import interp1d

def IsInside_CI(CI,mu,funt):
    # tomo n datos para calcular mi estadístico t y luego usar el cinturón
    # para hallar el CI
    xs = np.random.exponential(mu, n)
    tp = funt(xs)
    print(tp)
    CI = getCI(CI_t, taus, tp)
    return (CI[0] <= mu and mu <= CI[1])

def getCI(CI_t, taus, tp):
    # interpolo linealmente para hallar el CI
    f_l = interp1d(CI_t[:,1], taus); f_u = interp1d(CI_t[:,0], taus)
    tau_l = f_l(tp); tau_u = f_u(tp)
    return [tau_l, tau_u]

CI = getCI(CI_t, taus, tp)

plt.figure()    
plt.plot(CI_t[:,0], taus,'.b',CI_t[:,1], taus, '.b')
plt.vlines(tp, min(taus), max(taus), 'r')
plt.hlines(CI[0], np.min(CI_t), np.max(CI_t), 'g', alpha=0.8)
plt.hlines(CI[1], np.min(CI_t), np.max(CI_t), 'g', alpha=0.8)
plt.xlabel(r'$t$')
plt.ylabel(r'$\tau$')
plt.grid()
#%% ####IN PROCESS##### cobertura vs tau para las partes A y B ###############

# tenemos que alejanos de los límites de la interpolación para poder usarla
mu_min = 0.4
mu_max = 8.0
nscan_points= 1000
mus = np.linspace(mu_min, mu_max,nscan_points)

probs = cobertura_vs_par(mus,IsInside=IsInside_CI,funt=funt) 

plt.figure()    
plt.plot(mus,probs,'-r')
plt.hlines(CL, min(mus), max(mus), ls='dashed', lw = 0.8, alpha=0.8)
plt.xlabel(r'$\mu$')
plt.ylabel(r'$CL$')
plt.grid()
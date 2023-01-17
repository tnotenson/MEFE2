#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:46:16 2021

@author: tomasnotenson
"""
import numpy as np
import matplotlib.pyplot as plt
from time import time
from scipy.stats import pearsonr

plt.rcParams.update({
"text.usetex": True,
"font.family": "sans-serif",
"font.sans-serif": ["Helvetica"], "font.size": 20})

# Parte A
# 1. Agregar calculo de bias al script bootstrap_vs_toy que usamos ayer
# 2. Estudiar cómo evoluciona la estimación de la incerteza de un estimador a 
# elección conforme aumenta el número de réplicas para tres tamaños de 
# muestra diferente. Comparar en un mismo grafico.

# Parte B 
# Calcular la correlación angular en un experimento con fotones entrelazados
# usando los datos disponibles en la solapa Material Adicional, y agregar una barra
# de error por bootstrap para cada delta theta.

# A)
def bootstrap_poisson(estadistico, sample_size, replicas, CL=0.6827, mu=50, error='frecuentista'):
    # Arrays
    sample = np.random.poisson(mu, sample_size)
    bootstrap_sample = np.zeros(replicas)

    # Bootstrap
    for i in range(replicas):
        bootstrap_sample[i] = estadistico(np.random.choice(sample, size=sample_size, replace=True))

    # Bias
    bootstrap_mean = np.mean(bootstrap_sample)
    sample_statistic = estadistico(sample)
    bootstrap_bias = abs(bootstrap_mean-sample_statistic)

    # Sigma
    bootstrap_sigma = errores(bootstrap_sample, CL, tipo = error)
    
    return bootstrap_sigma, bootstrap_bias

def errores(sample, CL, tipo='std'):
    if tipo == 'std':
        sigma = np.std(sample, ddof=1)
        
    
    elif tipo == 'frecuentista':
        alpha = (1-CL)/2
        q = np.quantile(sample, [alpha, 1-alpha])
        sigma = (q[1]-q[0])/2
    
    return sigma
    
#%%
# Variables a probar
muestras = [20, 200, 2000] # número de muestras a probar
n_replicas = 200 # número de réplicas
replicas = 2*np.logspace(1, 4, n_replicas, dtype=int) # pido con los valores sean enteros
error = 'frecuentista' #tipo de error

# Estadístico
estadistico = np.mean
s_est = 'mean'

# estadistico = np.median
# s_est = 'median'

# estadistico = np.max
# s_est = 'max'

# estadistico = np.min
# s_est = 'min'

# estadistico = (lambda x: np.std(x, ddof=1))
# s_est = 'std'
                      
# CL = 0.6827; alpha = 1 - CL
# estadistico = (lambda x: np.quantile(x, q=[alpha]))
# s_est = 'quantile two tails CL = 68.27%'

# Itero sobre cada muestra
resultados_sigma, resultados_bias = {}, {}

start_time = time()

for sample_size in muestras:
    print(sample_size)
    sigma, bias = np.zeros(n_replicas), np.zeros(n_replicas)
    for i, r in enumerate(replicas):
        s, b = bootstrap_poisson(estadistico, sample_size, r)
        sigma[i], bias[i] = s, b

    resultados_sigma[sample_size] = sigma
    resultados_bias[sample_size] = bias

end_time = time()

print(f"Tardó: {end_time-start_time:.2f}s")

#%%
colors = ['r.', 'b.', 'g.']

fig = plt.figure(figsize=(12,8))
plt.title(f'Estadistico: {s_est}')
for k, sample_size in enumerate(muestras):
    plt.plot(replicas, resultados_sigma[sample_size], colors[k], 
            label=f'Tamaño = {sample_size}')
plt.xscale('log')
plt.xlabel('Replicas')
plt.ylabel(r'$\sigma$')
plt.legend(loc='best')
plt.grid()
plt.savefig('sigma_'+s_est+'.png')

#%% B)
fotones = np.array([49,57,50,54,65,45,55,52,57,61,47,53,44,57,57,53,39,58,47,61,50,60,63,44,57,46,59,61,46,57,44, 55,54,51,43,58])

def traslacion(fotones, n):
    return np.concatenate([fotones[-n:],fotones[:-n]])

def corr(a, b):
    c, _ = pearsonr(a,b) # me quedo con la correlación y tiro el p valor
    return c

def corr_mean(fotones):
    N = len(fotones)
    correlacion = np.zeros(N)
    for n in range(N):
        correlacion[n] = corr(fotones, traslacion(fotones, n)) 
    return correlacion

def boostrap_fotones(fotones, n, n_replicas=100, CL = 0.68, error='frecuentista'):
    # Cantidad de datos
    sample_size = len(fotones)

    # Calculo el array de conteos trasladado n
    traslacion_fotones  = traslacion(fotones, n)

    # Array de correlaciones
    correlacion = np.zeros(n_replicas)

    for i in range(n_replicas):
        # Selección random de sample_size muestras
        idx = np.random.choice(np.arange(sample_size), size=sample_size, replace=True)

        # Calculo correlación y almaceno
        correlacion[i] = corr(fotones[idx], traslacion_fotones[idx])

    bias = np.mean(correlacion)-corr(fotones, traslacion_fotones)

    # Sigma
    sigma = errores(correlacion, CL, tipo=error)

    return sigma, bias, np.mean(correlacion)

def bootstrap_correlacion(fotones, n_replicas, CL = 0.68, error='frecuentista'):
    # Variables
    N = len(fotones)
    correlacion_mean = np.zeros(N)
    correlacion_std = np.zeros(N)
    correlacion_bootstrap = np.zeros(N)
    correlacion_bias = np.zeros(N)
    
    # Lleno los arrays para cada traslación n
    for n in range(N):
        correlacion_mean[n] = corr(fotones, traslacion(fotones,n))
        sigma, bias, c_bootstrap = boostrap_fotones(fotones, n, n_replicas=n_replicas)
        correlacion_std[n], correlacion_bootstrap[n], correlacion_bias[n] = sigma, c_bootstrap, bias

    return correlacion_std, correlacion_bootstrap, correlacion_bias
#%%
n_replicas = 10000

std, mean, bias = bootstrap_correlacion(n_replicas)
correlacion = corr_mean(fotones)

#%%
x = np.arange(36)

plt.figure(figsize=(12,8))
plt.plot(x, correlacion, '^k')
plt.errorbar(x, correlacion, yerr=std, fmt='k--', ecolor='r')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$\rho$')
plt.grid()
plt.savefig("correlacion.png")

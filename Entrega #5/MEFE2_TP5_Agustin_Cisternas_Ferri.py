#%%
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
import pickle
from time import time
from scipy.stats import pearsonr

from numpy.random.mtrand import sample
from decimal import *
import array

# plt.rcParams.update({
# "text.usetex": True,
# "font.family": "sans-serif",
# "font.sans-serif": ["Helvetica"]})

##################################################################################################################################################################

"""
Parte A
1. Agregar calculo de bias al script bootstrap_vs_toy que usamos ayer
2. Estudiar cómo evoluciona la estimación de la incerteza de un estimador a 
elección conforme aumenta el número de réplicas para tres tamaños de 
muestra diferente. Comparar en un mismo grafico.
"""
#%%
def bootstrap_poisson(estadistico, sample_size, replicas, CL=0.68, mu=50, error='frecuentista'):
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
    if error == 'std':
        bootstrap_sigma = np.std(bootstrap_sample, ddof=1)

    elif error == 'frecuentista':
        alpha = (1-CL)/2
        q = np.quantile(bootstrap_sample, [alpha, 1-alpha])
        bootstrap_sigma = (q[1]-q[0])/2
    
    return bootstrap_sigma, bootstrap_bias

# Variables a probar
muestras = [20, 200, 2000]
n_replicas = 200
replicas = 2*np.logspace(1, 4, n_replicas, dtype=int)
error = 'frecuentista'

# Estadístico
CL = 0.6827
alpha = 1 - CL

estadistico = (lambda x: np.quantile(x, q=[alpha]))
s_est = 'quantileCL68.3'

# Itero sobre cada muestra
resultados_sigma, resultados_bias = {}, {}

#%%

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

colors = ['r.', 'g.', 'b.']

fig = plt.figure(figsize=(16,8))
for k, sample_size in enumerate(muestras):
    plt.plot(replicas, resultados_sigma[sample_size], colors[k], 
            label=f'Tamaño = {sample_size}')
plt.xscale('log')
plt.xlabel('Replicas')
plt.ylabel('Sigma')
plt.grid()
plt.savefig('sigma_'+s_est+'.png')


##################################################################################################################################################################
#%%
# fotones = np.array([49,57,50,54,65,45,55,52,57,61,47,53,
# 44,57,57,53,39,58,47,61,50,60,63,44,57,46,59,61,46,57,44,
# 55,54,51,43,58])

# def traslacion(fotones, n):
#     return np.concatenate([fotones[-n:],fotones[:-n]])

# def corr_mean(fotones):
#     N = len(fotones)
#     correlacion = np.zeros(N)
#     for n in range(N):
#         correlacion[n], _ = pearsonr(fotones, traslacion(fotones, n))
#     return correlacion

# def corr(a, b):
#     c, _ = pearsonr(a,b)
#     return c

# # %%

# def boostrap_fotones(fotones, n, n_replicas=100, 
#                     CL = 0.68, error='frecuentista'):
#     # Cantidad de datos
#     N = len(fotones)
#     sample_size = N

#     # Calculo el array de conteos trasladado n
#     traslacion_fotones  = traslacion(fotones, n)

#     # Array de correlaciones
#     correlacion = np.zeros(n_replicas)

#     for i in range(n_replicas):
#         # Selección random de sample_size muestras
#         idx = np.random.choice(np.arange(N), size=sample_size)

#         # Calculo correlación y almaceno
#         correlacion[i] = corr(fotones[idx], traslacion_fotones[idx])

#     bias = np.mean(correlacion)-corr(fotones, traslacion_fotones)

#     # Sigma
#     if error == 'std':
#         sigma = np.std(correlacion, ddof=1)

#     elif error == 'frecuentista':
#         alpha = (1-CL)/2
#         q = np.quantile(correlacion, [alpha, 1-alpha])
#         sigma = (q[1]-q[0])/2

#     return sigma, bias, np.mean(correlacion)

# # %%
# def bootstrap_correlacion(n_replicas):
#     # Variables
#     N = len(fotones)
#     correlacion_mean = np.zeros(N)
#     correlacion_std = np.zeros(N)
#     correlacion_bootstrap = np.zeros(N)
#     correlacion_bias = np.zeros(N)

#     for n in range(N):
#         correlacion_mean[n] = corr(fotones, traslacion(fotones,n))
#         sigma, bias, c_bootstrap = boostrap_fotones(fotones, n, n_replicas=n_replicas)
#         correlacion_std[n], correlacion_bootstrap[n], correlacion_bias[n] = sigma, c_bootstrap, bias

#     return correlacion_std, correlacion_bootstrap, correlacion_bias
# #%%
# n_replicas = 10000

# std, mean, bias = bootstrap_correlacion(n_replicas)
# correlacion = corr_mean(fotones)

# #%%

# x = np.arange(36)

# plt.figure(figsize=(16,8))
# plt.errorbar(x, correlacion, yerr=std, fmt='k-', ecolor='r')
# plt.xlabel(r'$\theta$')
# plt.ylabel(r'$\rho$')
# plt.grid()
# plt.savefig("correlacion.png")

# %%

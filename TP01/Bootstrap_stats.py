#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 15:45:40 2017

@author: Charles
"""

import numpy as np
import datetime
import matplotlib.pyplot as plt

### Input data for pole of plane 1
azi_pole_plane1 = np.array([311,310,320,302,300,298,300,300,312,294,302,303,319,320,305,310,310,302,313])
dip_pole_plane1 = np.array([03,20,25,15,10,12,26,16,29,29,16,23,23,18,23,25,35,16,22])
pole_plane1 = np.vstack((azi_pole_plane1, dip_pole_plane1)).T
### Input data for pole of plane 2
azi_pole_plane2 = np.array([190,144,128,82,150,162,146,160,102,150,132,160,114,135,150,143,140,135,179])
dip_pole_plane2 = np.array([84,69,64,71,78,74,62,70,57,56,74,62,64,72,65,64,55,74,60])
pole_plane2 = np.vstack((azi_pole_plane2, dip_pole_plane2)).T
### Input data for P axis
azi_P_axis = np.array([136,134,136,112,125,127,129,130,121,128,124,135,131,139,133,134,134,125,148])
plu_P_axis = np.array([42,25,20,29,34,32,19,28,15,14,29,20,21,27,21,19,10,29,20])
P_axis = np.vstack((azi_P_axis, plu_P_axis)).T

### Nb of observations and iterations for bootstrap MonteCarlo statistics                  
nb_data = len(P_axis)
iterations = 10000

### Empty arrays to store results
means_pole_plane1 = np.empty((iterations,2))
means_pole_plane2 = np.empty((iterations,2))
means_P_axis = np.empty((iterations,2))

### Bootstrap statistics
for i in range(iterations):
      # Pick 9 random data from the data set
      pick = np.random.randint(0, nb_data, size=nb_data/2)
      # Compute mean of randomly picked data for plane1, plane2, P_axis
      means_pole_plane1[i] = np.mean(pole_plane1[pick], axis=0)
      means_pole_plane2[i] = np.mean(pole_plane2[pick], axis=0)
      means_P_axis[i] = np.mean(P_axis[pick], axis=0)
      
### Plots
fig, ax = plt.subplots(2,1, figsize=(8,6))
ax[0].hist(pole_plane1[:,0], label="Azimuth ($^o$)")
ax[1].hist(pole_plane1[:,1], label="Dip ($^o$)")
ax[0].legend(), ax[1].legend()
ax[0].set_ylabel("Occurence"), ax[1].set_ylabel("Occurence")
plt.suptitle("Data distribution for pole of plane 1", fontsize=14)
fig.savefig("Data_plane1")
      
fig, ax = plt.subplots(2,1, figsize=(8,6))
ax[0].hist(means_pole_plane1[:,0], bins=20, label="Azimuth ($^o$)\n$\mu$: %.1f\n$\sigma$: %.1f" %(np.mean(means_pole_plane1[:,0]),np.std(means_pole_plane1[:,0])))
ax[1].hist(means_pole_plane1[:,1], bins=20, label="Dip ($^o$)\n$\mu$: %.1f\n$\sigma$: %.1f" %(np.mean(means_pole_plane1[:,1]),np.std(means_pole_plane1[:,1])))
ax[0].legend(), ax[1].legend()
ax[0].set_ylabel("Occurence"), ax[1].set_ylabel("Occurence")
plt.suptitle("Bootstrap distribution for pole of plane 1", fontsize=14)
fig.savefig("Bootstrap_plane1")

fig, ax = plt.subplots(2,1, figsize=(8,6))
ax[0].hist(pole_plane2[:,0], label="Azimuth ($^o$)")
ax[1].hist(pole_plane2[:,1], label="Dip ($^o$)")
ax[0].legend(), ax[1].legend()
ax[0].set_ylabel("Occurence"), ax[1].set_ylabel("Occurence")
plt.suptitle("Data distribution for pole of plane 2", fontsize=14)
fig.savefig("Data_plane2")
      
fig, ax = plt.subplots(2,1, figsize=(8,6))
ax[0].hist(means_pole_plane2[:,0], bins=20, label="Azimuth ($^o$)\n$\mu$: %.1f\n$\sigma$: %.1f" %(np.mean(means_pole_plane2[:,0]),np.std(means_pole_plane2[:,0])))
ax[1].hist(means_pole_plane2[:,1], bins=20, label="Dip ($^o$)\n$\mu$: %.1f\n$\sigma$: %.1f" %(np.mean(means_pole_plane2[:,1]),np.std(means_pole_plane2[:,1])))
ax[0].legend(), ax[1].legend()
ax[0].set_ylabel("Occurence"), ax[1].set_ylabel("Occurence")
plt.suptitle("Bootstrap distribution for pole of plane 2", fontsize=14)
fig.savefig("Bootstrap_plane2")

fig, ax = plt.subplots(2,1, figsize=(8,6))
ax[0].hist(means_P_axis[:,0], label="Azimuth ($^o$)")
ax[1].hist(means_P_axis[:,1], label="Dip ($^o$)")
ax[0].legend(), ax[1].legend()
ax[0].set_ylabel("Occurence"), ax[1].set_ylabel("Occurence")
plt.suptitle("Data distribution for P axis", fontsize=14)
fig.savefig("Data_Paxis")
      
fig, ax = plt.subplots(2,1, figsize=(8,6))
ax[0].hist(means_P_axis[:,0], bins=20, label="Azimuth ($^o$)\n$\mu$: %.1f\n$\sigma$: %.1f" %(np.mean(means_P_axis[:,0]),np.std(means_P_axis[:,0])))
ax[1].hist(means_P_axis[:,1], bins=20, label="Plunge ($^o$)\n$\mu$: %.1f\n$\sigma$: %.1f" %(np.mean(means_P_axis[:,1]),np.std(means_P_axis[:,1])))
ax[0].legend(), ax[1].legend()
ax[0].set_ylabel("Occurence"), ax[1].set_ylabel("Occurence")
plt.suptitle("Bootstrap distribution for P axis", fontsize=14)
fig.savefig("Bootstrap_Paxis")
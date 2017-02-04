# -*- coding: utf-8 -*-
"""
Created on Sat Feb 04 10:23:52 2017

@author: Charles
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

delta_T0_diurne = 20 # C
theta_diurne = 24*60*60 # s
delta_T0_annuel = 10 # C
theta_annuel = 365*24*60*60 # s
kappa = (1e-4)*4e-3 # m^2 s^-1
z = np.linspace(0,2.5,100) # 0 to 2.5 meters, 100 points
t_diurne = np.linspace(0,24*60*60,100) # En secondes
t_annuel = np.linspace(0,365*24*60*60,100) # En secondes
             
def var_temp_max(delta_T0, z, kappa, theta):
    delta_T = delta_T0 * np.exp(-z/(np.sqrt(kappa*theta/np.pi)))
    return delta_T # C

def var_temp_temp(delta_T0, z, kappa, theta, t):
    delta_T = delta_T0 * np.exp(-z/(np.sqrt(kappa*theta/np.pi))) * np.sin((2*np.pi*t/theta) - (z/(np.sqrt(kappa*theta/np.pi))))
    return delta_T # C

# Variation dirune de 20 degrés
var_diurne = var_temp_max(delta_T0_diurne, z, kappa, theta_diurne)
var40cm_diurne = var_temp_max(delta_T0_diurne, 0.4, kappa, theta_diurne)
var2m_diurne = var_temp_max(delta_T0_diurne, 2, kappa, theta_diurne)

# Variation annuel de 10 degrés
var_annuel = var_temp_max(delta_T0_annuel, z, kappa, theta_annuel)
var40cm_annuel = var_temp_max(delta_T0_annuel, 0.4, kappa, theta_annuel)
var2m_annuel = var_temp_max(delta_T0_annuel, 2, kappa, theta_annuel)

# Temps requis pour que la variation maximale soit atteinte à 40 cm
var_temps_40cm_diurne = var_temp_temp(delta_T0_diurne, 0.40, kappa, theta_diurne, t_diurne)
var_temps_40cm_annuel = var_temp_temp(delta_T0_annuel, 0.40, kappa, theta_annuel, t_annuel)

temps_requis_diurne = t_diurne[np.argmin(np.abs(var_temps_40cm_diurne - var40cm_diurne))] /(60*60) # En heure
temps_requis_annuel = t_annuel[np.argmin(np.abs(var_temps_40cm_annuel - var40cm_annuel))] /(60*60*24) # En jours

print var40cm_diurne, var2m_diurne
print var40cm_annuel, var2m_annuel
print temps_requis_diurne
print temps_requis_annuel

# Figure
fig1, ax = plt.subplots(1,2,figsize=(8,3))
ax[0].semilogx(var_diurne,z,'-b')
ax[0].plot(var40cm_diurne,0.4,'ok')
ax[0].plot(var2m_diurne,2,'ok')
ax[1].plot(var_annuel,z,'-g')
ax[1].plot(var40cm_annuel,0.4,'ok')
ax[1].plot(var2m_annuel,2,'ok')
ax[0].set_xlabel("Variation temperature maximale ($^o$C)")
ax[1].set_xlabel("Variation temperature maximale ($^o$C)")
ax[0].set_ylabel("Profondeur (m)")
ax[0].invert_yaxis(), ax[1].invert_yaxis()
plt.legend()
ax[0].set_title("Variation diurne de 20$^o$C")
ax[1].set_title("Variation annuelle de 10$^o$C")
plt.tight_layout()

fig2, ax = plt.subplots(1,2,figsize=(8,3))
ax[0].plot(t_diurne/(60*60),var_temps_40cm_diurne,'-b')
ax[0].plot(temps_requis_diurne,var40cm_diurne,'ok')
ax[1].plot(t_annuel/(60*60*24),var_temps_40cm_annuel,'-g')
ax[1].plot(temps_requis_annuel,var40cm_annuel,'ok')
ax[0].set_xlabel("Temps (heures)")
ax[1].set_xlabel("Temps (jours)")
ax[0].set_ylabel("Variation temperature ($^o$C)")
plt.legend()
ax[0].set_title(u"Variation diurne de 20$^o$C à 40 cm")
ax[1].set_title(u"Variation annuelle de 10$^o$C à 40 cm")
plt.tight_layout()


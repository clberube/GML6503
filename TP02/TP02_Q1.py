# -*- coding: utf-8 -*-
"""
Created on Sat Feb 04 09:09:57 2017

@author: Charles
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

A_0 = 12.6e-7 # W m^-3 
dT_0 = 20*(1e-3) # C m^1
k = 1.68 # W m^-1 C^-1 
T_0 = 0 # C
Z_eval = 40*(1e3) # m
z = np.linspace(0,50,100)*1000 # z from 0 to 80, 100 points total
    
def gradient_thermique(A_0, z, k, dT_0, D=20*(1e3), a_cte=False):
    if a_cte:
        return -A_0*z/k  + dT_0 # Returns in C m^-1
    else:
        return np.exp(-z/D)*A_0*D/k  + dT_0 - A_0*D/k # Returns in C m^-1

def temperature(A_0, z, k, dT_0, D=20*(1e3), a_cte=False):
    if a_cte:
        return (-0.5*A_0/k)*z**2  + dT_0*z # Returns in C
    else:
        return -1*np.exp(-z/D)*A_0*(D**2)/k  + z*(dT_0 - A_0*D/k) + A_0*(D**2)/k # Returns in C m^-1

temp_a_cte = temperature(A_0,Z_eval,k,dT_0,a_cte=True) # In C
grad_a_cte = 1000*gradient_thermique(A_0,Z_eval,k,dT_0,a_cte=True) # In C km^-1
print temp_a_cte, grad_a_cte

temp_a = temperature(A_0,Z_eval,k,dT_0) # In C
grad_a = 1000*gradient_thermique(A_0,Z_eval,k,dT_0) # In C km^-1
print temp_a, grad_a

# Figures
fig1, ax1 = plt.subplots(1,2,figsize=(8,3))
ax1[0].plot(1000*gradient_thermique(A_0,z,k,dT_0,a_cte=True),z/1000,'-b')
ax1[0].plot(grad_a_cte,Z_eval/1000,'ok')
ax1[1].plot(temperature(A_0,z,k,dT_0,a_cte=True),z/1000,'-g')
ax1[1].plot(temp_a_cte,Z_eval/1000,'ok')
ax1[0].set_xlabel("Gradient thermique ($^o$C/km)")
ax1[1].set_xlabel("Temperature ($^o$C)")
ax1[0].set_ylabel("Profondeur (km)")
ax1[0].invert_yaxis(), ax1[1].invert_yaxis()
plt.legend()
ax1[0].set_title("$A = A_0 = 12.6 x 10 ^{-7} \/ W/m^3$")
ax1[1].set_title("$A = A_0 = 12.6 x 10 ^{-7} \/ W/m^3$")
plt.tight_layout()

fig2, ax2 = plt.subplots(1,2,figsize=(8,3))
ax2[0].plot(1000*gradient_thermique(A_0, z, k, dT_0),z/1000,'-b')
ax2[0].plot(grad_a,Z_eval/1000,'ok')
ax2[1].plot(temperature(A_0,z,k,dT_0),z/1000,'-g')
ax2[1].plot(temp_a,Z_eval/1000,'ok')
ax2[0].set_xlabel("Gradient thermique ($^o$C/km)")
ax2[1].set_xlabel("Temperature ($^o$C)")
ax2[0].set_ylabel("Profondeur (km)")
ax2[0].invert_yaxis(), ax2[1].invert_yaxis()
ax2[0].set_title("$A = A_0 e^{-z/D} \/ W/m^3$")
ax2[1].set_title("$A = A_0 e^{-z/D} \/ W/m^3$")
plt.legend()
plt.tight_layout()
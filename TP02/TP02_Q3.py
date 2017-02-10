#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 16:25:12 2017

@author: Charles
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

A = {"Grani": 1.0*1000, # En W km^-3
     "Metas": 0.8*1000,
     "Granu": 0.2*1000,
     "Perid": 0.0*1000,}

k = {"Grani": 3.5*1000, # En W km^-1
     "Metas": 4.0*1000,
     "Granu": 3.0*1000,
     "Perid": 2.5*1000,}

rho = {"Grani": 2.60*1e3,
       "Metas": 2.65*1e3,
       "Granu": 2.90*1e3,
       "Perid": 3.30*1e3,}

B = {"Grani": 1.80e-3,
     "Metas": 1.26e-3,
     "Granu": 3.10e-5,
     "Perid": 1.61e+2,}

n = {"Grani": 2.0,
     "Metas": 2.4,
     "Granu": 3.5,
     "Perid": 3.4,}

Q = {"Grani": 156*1000,
     "Metas": 180*1000,
     "Granu": 236*1000,
     "Perid": 474*1000,}

Abitibi = {"Q": 38,
           "z1": 10,
           "z2": 25,
           "z3": 35,
           "z4": 150,}

Grenville = {"Q": 40,
             "z1": 10,
             "z2": 30,
             "z3": 40,
             "z4": 150,}

Appalaches = {"Q": 57,
              "z1": 15,
              "z2": 30,
              "z3": 45,
              "z4": 150,}

Cordillera = {"Q": 80,
              "z1": 20,
              "z2": 30,
              "z3": 40,
              "z4": 150,}

def frottement(rho,z,g=9.81,mu=0.6):
    sigma = mu*rho*g*z
    return sigma

def fluage(B, n, Q, T, epsilon=1e-14, R=8.314):
    sigma = ((epsilon/B)**(1/n))*np.exp(Q/(n*R*T))
    return sigma

def C1(Q0,k1):
    return Q0/k1

def C3(Q0,A1,A2,k1,k2,z1):
    return C1(Q0,k1) + z1*(A2/k2 - A1/k1)
    
def C4(A1,A2,k1,k2,z1):
    return (z1**2) * (A1/k1 - A2/k2)/2
    
def C5(Q0,A1,A2,A3,k1,k2,k3,z1,z2):
    return z2*(A3/k3 - A2/k2) + z1*(A2/k2 - A1/k1) + C1(Q0,k1)

def C6(A1,A2,A3,k1,k2,k3,z1,z2):
    return ((z2**2)*(A2/k2 - A3/k3) + (z1**2)*(A1/k1 - A2/k2))/2

def C7(Q0,A1,A2,A3,A4,k1,k2,k3,k4,z1,z2,z3):
    return z3*(A4/k4 - A3/k3) + z2*(A3/k3 - A2/k2) + z1*(A2/k2 - A1/k1) + C1(Q0,k1)
    
def C8(A1,A2,A3,A4,k1,k2,k3,k4,z1,z2,z3):
    return ((z3**2)*(A3/k3 - A4/k4) + (z2**2)*(A2/k2 - A3/k3) + (z1**2)*(A1/k1 - A2/k2))/2

def temperature(region,A,k):
    Q0 = region["Q"]*1000 # En W km^-2
    z1 = np.arange(0,region["z1"]+1)
    z2 = np.arange(region["z1"],region["z2"]+1)
    z3 = np.arange(region["z2"],region["z3"]+1)
    z4 = np.arange(region["z3"],region["z4"]+1)

    frot_grani=frottement(rho["Grani"], z1*1000)*(1e-6) # In MPa
    frot_metas=frottement(rho["Metas"], z2*1000)*(1e-6) # In MPa
    frot_granu=frottement(rho["Granu"], z3*1000)*(1e-6) # In MPa
    frot_perid=frottement(rho["Perid"], z4*1000)*(1e-6) # In MPa
    
    T_grani = -A["Grani"]*(z1**2)/(2*k["Grani"]) + z1*C1(Q0,k["Grani"])
    dT_grani = -A["Grani"]*(z1)/(k["Grani"]) + C1(Q0,k["Grani"])

    T_metas = -A["Metas"]*(z2**2)/(2*k["Metas"]) + z2*C3(Q0,A["Grani"],A["Metas"],k["Grani"],k["Metas"],region["z1"]) + C4(A["Grani"],A["Metas"],k["Grani"],k["Metas"],region["z1"])
    dT_metas = -A["Metas"]*(z2)/(k["Metas"]) + C3(Q0,A["Grani"],A["Metas"],k["Grani"],k["Metas"],region["z1"])

    
    T_granu = -A["Granu"]*(z3**2)/(2*k["Granu"]) + z3*C5(Q0,A["Grani"],A["Metas"],A["Granu"],k["Grani"],k["Metas"],k["Granu"],region["z1"],region["z2"]) + C6(A["Grani"],A["Metas"],A["Granu"],k["Grani"],k["Metas"],k["Granu"],region["z1"],region["z2"])
    dT_granu = -A["Granu"]*(z3)/(k["Granu"]) + C5(Q0,A["Grani"],A["Metas"],A["Granu"],k["Grani"],k["Metas"],k["Granu"],region["z1"],region["z2"])

    
    
    T_perid = -A["Perid"]*(z4**2)/(2*k["Perid"]) + z4*C7(Q0,A["Grani"],A["Metas"],A["Granu"],A["Perid"],k["Grani"],k["Metas"],k["Granu"],k["Perid"],region["z1"],region["z2"],region["z3"]) + C8(A["Grani"],A["Metas"],A["Granu"],A["Perid"],k["Grani"],k["Metas"],k["Granu"],k["Perid"],region["z1"],region["z2"],region["z3"])
    dT_perid = -A["Perid"]*(z4)/(k["Perid"]) + C7(Q0,A["Grani"],A["Metas"],A["Granu"],A["Perid"],k["Grani"],k["Metas"],k["Granu"],k["Perid"],region["z1"],region["z2"],region["z3"])
    
    
    
    flua_grani=fluage(B["Grani"], n["Grani"], Q["Grani"], T_grani+273)
    flua_metas=fluage(B["Metas"], n["Metas"], Q["Metas"], T_metas+273)
    flua_granu=fluage(B["Granu"], n["Granu"], Q["Granu"], T_granu+273)
    flua_perid=fluage(B["Perid"], n["Perid"], Q["Perid"], T_perid+273)
    
    return np.hstack((z1, z2, z3, z4)), np.hstack((T_grani, T_metas, T_granu, T_perid)), np.hstack((dT_grani, dT_metas, dT_granu, dT_perid)), np.hstack((frot_grani, frot_metas, frot_granu, frot_perid)), np.hstack((flua_grani, flua_metas, flua_granu, flua_perid))

plt.rcParams['savefig.dpi'] = 72
plt.rcParams['figure.dpi'] = 72
fig1, ax1 = plt.subplots(2,2,figsize=(10,8))
fig2, ax2 = plt.subplots(2,2,figsize=(10,8))
fig3, ax3 = plt.subplots(2,2,figsize=(10,8))

z, T, dT, fr, fl = temperature(Abitibi,A,k)
ax1[0,0].plot(T,z,"k-",linewidth=2)
ax1[0,0].set_xlabel("Temperature ($^o$C)")
ax1[0,0].set_ylabel("Profondeur (km)")
ax1[0,0].set_title("Abitibi")
ax1[0,0].invert_yaxis() 
ax3[0,0].plot(dT,z, "r-",linewidth=2)
ax3[0,0].set_xlabel("Gradient thermique ($^o$C / km)")
ax3[0,0].set_ylabel("Profondeur (km)")
ax3[0,0].set_title("Abitibi")
ax3[0,0].invert_yaxis() 
fr2 = np.ma.masked_greater(fr, fl)
fl2 = np.ma.masked_greater(fl, fr)
ax2[0,0].semilogx(fr2,z,'b-',label="Frottement",linewidth=2)
ax2[0,0].semilogx(fl2,z,'g-',label="Fluage",linewidth=2)
ax2[0,0].semilogx(fr,z,'b-',alpha=0.2)
ax2[0,0].semilogx(fl,z,'g-',alpha=0.2)
ax2[0,0].legend(loc=4)
ax1[0,0].axhline(y=Abitibi["z1"],color='grey',zorder=-1)
ax1[0,0].axhline(y=Abitibi["z2"],color='grey',zorder=-1)
ax1[0,0].axhline(y=Abitibi["z3"],color='grey',zorder=-1)
ax1[0,0].axhline(y=Abitibi["z4"],color='grey',zorder=-1)
ax2[0,0].axhline(y=Abitibi["z1"],color='grey',zorder=-1)
ax2[0,0].axhline(y=Abitibi["z2"],color='grey',zorder=-1)
ax2[0,0].axhline(y=Abitibi["z3"],color='grey',zorder=-1)
ax2[0,0].axhline(y=Abitibi["z4"],color='grey',zorder=-1)
ax3[0,0].axhline(y=Abitibi["z1"],color='grey',zorder=-1)
ax3[0,0].axhline(y=Abitibi["z2"],color='grey',zorder=-1)
ax3[0,0].axhline(y=Abitibi["z3"],color='grey',zorder=-1)
ax3[0,0].axhline(y=Abitibi["z4"],color='grey',zorder=-1)
ax2[0,0].set_xlabel("$\sigma$ (MPa)")
ax2[0,0].set_ylabel("Profondeur (km)")
ax2[0,0].set_title("Abitibi")
ax2[0,0].invert_yaxis() 


z, T, dT, fr, fl = temperature(Grenville,A,k)
ax1[0,1].plot(T,z,"k-",linewidth=2)
ax1[0,1].set_xlabel("Temperature ($^o$C)")
ax1[0,1].set_ylabel("Profondeur (km)")
ax1[0,1].set_title("Grenville")
ax1[0,1].invert_yaxis() 
ax3[0,1].plot(dT,z, "r-",linewidth=2)
ax3[0,1].set_xlabel("Gradient thermique ($^o$C / km)")
ax3[0,1].set_ylabel("Profondeur (km)")
ax3[0,1].set_title("Grenville")
ax3[0,1].invert_yaxis() 
fr2 = np.ma.masked_greater(fr, fl)
fl2 = np.ma.masked_greater(fl, fr)
ax2[0,1].semilogx(fr2,z,'b-',label="Frottement",linewidth=2)
ax2[0,1].semilogx(fl2,z,'g-',label="Fluage",linewidth=2)
ax2[0,1].semilogx(fr,z,'b-',alpha=1)
ax2[0,1].semilogx(fl,z,'g-',alpha=1)
ax2[0,1].legend(loc=4)
ax1[0,1].axhline(y=Grenville["z1"],color='grey',zorder=-1)
ax1[0,1].axhline(y=Grenville["z2"],color='grey',zorder=-1)
ax1[0,1].axhline(y=Grenville["z3"],color='grey',zorder=-1)
ax1[0,1].axhline(y=Grenville["z4"],color='grey',zorder=-1)
ax2[0,1].axhline(y=Grenville["z1"],color='grey',zorder=-1)
ax2[0,1].axhline(y=Grenville["z2"],color='grey',zorder=-1)
ax2[0,1].axhline(y=Grenville["z3"],color='grey',zorder=-1)
ax2[0,1].axhline(y=Grenville["z4"],color='grey',zorder=-1)
ax3[0,1].axhline(y=Grenville["z1"],color='grey',zorder=-1)
ax3[0,1].axhline(y=Grenville["z2"],color='grey',zorder=-1)
ax3[0,1].axhline(y=Grenville["z3"],color='grey',zorder=-1)
ax3[0,1].axhline(y=Grenville["z4"],color='grey',zorder=-1)
ax2[0,1].set_xlabel("$\sigma$ (MPa)")
ax2[0,1].set_ylabel("Profondeur (km)")
ax2[0,1].set_title("Grenville")
ax2[0,1].invert_yaxis() 

z, T, dT, fr, fl = temperature(Appalaches,A,k)
ax1[1,0].plot(T,z,"k-",linewidth=2)
ax1[1,0].set_xlabel("Temperature ($^o$C)")
ax1[1,0].set_ylabel("Profondeur (km)")
ax1[1,0].set_title("Appalaches")
ax1[1,0].invert_yaxis() 
ax3[1,0].plot(dT,z, "r-",linewidth=2)
ax3[1,0].set_xlabel("Gradient thermique ($^o$C / km)")
ax3[1,0].set_ylabel("Profondeur (km)")
ax3[1,0].set_title("Appalaches")
ax3[1,0].invert_yaxis() 
fr2 = np.ma.masked_greater(fr, fl)
fl2 = np.ma.masked_greater(fl, fr)
ax2[1,0].semilogx(fr2,z,'b-',label="Frottement",linewidth=2)
ax2[1,0].semilogx(fl2,z,'g-',label="Fluage",linewidth=2)
ax2[1,0].semilogx(fr,z,'b-',alpha=0.2)
ax2[1,0].semilogx(fl,z,'g-',alpha=0.2)
ax2[1,0].legend(loc=4)
ax1[1,0].axhline(y=Appalaches["z1"],color='grey',zorder=-1)
ax1[1,0].axhline(y=Appalaches["z2"],color='grey',zorder=-1)
ax1[1,0].axhline(y=Appalaches["z3"],color='grey',zorder=-1)
ax1[1,0].axhline(y=Appalaches["z4"],color='grey',zorder=-1)
ax2[1,0].axhline(y=Appalaches["z1"],color='grey',zorder=-1)
ax2[1,0].axhline(y=Appalaches["z2"],color='grey',zorder=-1)
ax2[1,0].axhline(y=Appalaches["z3"],color='grey',zorder=-1)
ax2[1,0].axhline(y=Appalaches["z4"],color='grey',zorder=-1)
ax3[1,0].axhline(y=Appalaches["z1"],color='grey',zorder=-1)
ax3[1,0].axhline(y=Appalaches["z2"],color='grey',zorder=-1)
ax3[1,0].axhline(y=Appalaches["z3"],color='grey',zorder=-1)
ax3[1,0].axhline(y=Appalaches["z4"],color='grey',zorder=-1)
ax2[1,0].set_xlabel("$\sigma$ (MPa)")
ax2[1,0].set_ylabel("Profondeur (km)")
ax2[1,0].set_title("Appalaches")
ax2[1,0].invert_yaxis() 

z, T, dT, fr, fl = temperature(Cordillera,A,k)
ax1[1,1].plot(T,z,"k-",linewidth=2)
ax1[1,1].set_xlabel("Temperature ($^o$C)")
ax1[1,1].set_ylabel("Profondeur (km)")
ax1[1,1].set_title("Cordillera")
ax1[1,1].invert_yaxis()
ax3[1,1].plot(dT,z, "r-",linewidth=2)
ax3[1,1].set_xlabel("Gradient thermique ($^o$C / km)")
ax3[1,1].set_ylabel("Profondeur (km)")
ax3[1,1].set_title("Cordillera")
ax3[1,1].invert_yaxis() 
fr2 = np.ma.masked_greater(fr, fl)
fl2 = np.ma.masked_greater(fl, fr)
ax2[1,1].semilogx(fr2,z,'b-',label="Frottement",linewidth=2)
ax2[1,1].semilogx(fl2,z,'g-',label="Fluage",linewidth=2)
ax2[1,1].semilogx(fr,z,'b-',alpha=0.2)
ax2[1,1].semilogx(fl,z,'g-',alpha=0.2)
ax2[1,1].legend(loc=4)
ax2[1,1].axhline(y=Cordillera["z1"],color='grey',zorder=-1)
ax2[1,1].axhline(y=Cordillera["z2"],color='grey',zorder=-1)
ax2[1,1].axhline(y=Cordillera["z3"],color='grey',zorder=-1)
ax2[1,1].axhline(y=Cordillera["z4"],color='grey',zorder=-1)
ax1[1,1].axhline(y=Cordillera["z1"],color='grey',zorder=-1)
ax1[1,1].axhline(y=Cordillera["z2"],color='grey',zorder=-1)
ax1[1,1].axhline(y=Cordillera["z3"],color='grey',zorder=-1)
ax1[1,1].axhline(y=Cordillera["z4"],color='grey',zorder=-1)
ax3[1,1].axhline(y=Cordillera["z1"],color='grey',zorder=-1)
ax3[1,1].axhline(y=Cordillera["z2"],color='grey',zorder=-1)
ax3[1,1].axhline(y=Cordillera["z3"],color='grey',zorder=-1)
ax3[1,1].axhline(y=Cordillera["z4"],color='grey',zorder=-1)
ax2[1,1].set_xlabel("$\sigma$ (MPa)")
ax2[1,1].set_ylabel("Profondeur (km)")
ax2[1,1].set_title("Cordillera")
ax2[1,1].invert_yaxis()

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig1.savefig("Q3A.pdf")
fig2.savefig("Q3B.pdf")
fig3.savefig("Q3C.pdf")
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 23:50:41 2022

@author: utente
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from sklearn.metrics import r2_score
import pickle
from func_03_06_class import *


#%% using only tafel and semplify  everything
#%% setting conditions
T=353
COH_conc=4.212 #mol/dm^3
P=20 #oxygen pressure 20 atm
P_s=7 #bar system pressure bar
w=0.30
delta_l=0.0132+0.0350 #from thesis (separator and electrode thichness) expressed in [cm]

#%% reference conditions Miles experiments
# T_ref=353
# w_ref=0.35
# P_ref=20 #atm oxygen
#%%
a=[]
for i in np.linspace(0,30000,100):
    a.append(((i/30000)**(0.3)))
    
# for i in np.linspace(0,100000,100):
#     a.append(0.023*((i)**(0.3)))

plt.plot(np.linspace(0,30000,100),a)    
plt.ylim(0,1.0)

#%% load model class
model= AEC_model(T,COH_conc,P,P_s,w,delta_l)
#%%load reference conditions
# model.ref_conditions(T_ref,w_ref,P_ref)
# coh_text="$C_{OH,ref}$"
# co2_text="$C_{O2,ref}$"

# print(coh_text+f": {round(model.COH_ref,2)} $[mol/dm^3]$")
# print(co2_text+f": {round(model.CO2_ref,5)} $[mol/dm^3]$")


#%%


#%% multi pressure plots
# a=[]
# for i in np.arange(1,20,0.2):
#     a.append(model.v_rev(T,w,i))

# plt.figure(figsize=(15,7))
# plt.plot(np.arange(1,20,0.2),a)
# plt.xlabel("Operating pressure [bar]")
# plt.ylabel("$V_{rev}$ $[V]$")
# plt.title("Reversible voltage dependence")
# plt.grid()
#%% reaction distances
# plt.figure(235,figsize=(15,7))
# df_comsol=pd.read_csv("COMSOL.csv",index_col=0)
# df_exp=pd.read_csv("Experiment.csv",index_col=0)
# plt.plot([i/10 for i in df_comsol.index],df_comsol.iloc[:,0].tolist(),label="COMSOL model (thesis)",linewidth=4,c="#ff7f0e")
# plt.plot([i/10 for i in df_exp.index],df_exp.iloc[:,0].tolist(),"o",c='k',markeredgewidth=3, markersize=12,fillstyle="none", label="Experimental data",linewidth=3)  

# for i in np.linspace(0.0132+0.0350,0.0132+0.0350*2,5):
#     df_kibria_mix=model.kibria_mix(T,COH_conc,i,w,P_s,P)
#     plt.plot(df_kibria_mix.index,df_kibria_mix,label=f"Kibria Mix tafel with $\delta l$={round(i,3)} [cm]")

# plt.legend()
# plt.grid()
# plt.title("changing reaction distances [cm]")
    
    
#%% mixed kibria
df=model.kibria_mix(T,COH_conc,delta_l,w,P_s,P)

df_comsol=pd.read_csv("COMSOL.csv",index_col=0)
df_exp=pd.read_csv("Experiment.csv",index_col=0)

plt.figure(figsize=(18,9))

plt.plot([i/10 for i in df_comsol.index],df_comsol.iloc[:,0].tolist(),label="COMSOL model (thesis)",linewidth=4,c="#ff7f0e")
plt.plot([i/10 for i in df_exp.index],df_exp.iloc[:,0].tolist(),"o",c='k',markeredgewidth=3, markersize=12,fillstyle="none", label="Experimental data",linewidth=3)  
plt.plot(df.index,df,label="Kibria Mix Tafel")
plt.legend()
o="$P_{O2}$"
plt.title(f"Set-up: T={T} [K],COH={COH_conc} $[mol/dm^3]$,"+o+f"={P}[atm], w={100*w}[%],$P_s$={P_s} [bar]")
plt.tight_layout()
plt.grid()


#%%  moldel miles conditions
df_miles_ohm=model.Miles_mix_bub_ohm(T,COH_conc,delta_l,w,P_s,P)

df_miles=model.Miles_mix_bub_act(T,COH_conc,delta_l,w,P_s,P)

df_comsol=pd.read_csv("COMSOL.csv",index_col=0)
df_exp=pd.read_csv("Experiment.csv",index_col=0)

plt.figure(figsize=(18,9))

plt.plot([i/10 for i in df_comsol.index],df_comsol.iloc[:,0].tolist(),label="COMSOL model (thesis)",linewidth=4,c="r")
plt.plot([i/10 for i in df_exp.index],df_exp.iloc[:,0].tolist(),"o",c='k',markeredgewidth=3, markersize=12,fillstyle="none", label="Experimental data",linewidth=3)  
# plt.plot(df_no_bub.index,df_no_bub,label="Kibria Mix Tafel no bubbles activation")
plt.plot(df_miles_ohm.index,df_miles_ohm,label="Kibria Mix Tafel with bubbles ohmic")
plt.plot(df_miles.index,df_miles,label="Miles Mix Tafel without bubbles activation")


plt.legend()
o="$P_{O2}$"
plt.title(f"Set-up: T={T} [K],COH={COH_conc} $[mol/dm^3]$,"+o+f"={P}[atm], w={100*w}[%],$P_s$={P_s} [bar]")
plt.tight_layout()
plt.grid()

#%% bubble overvoltage
df_no_bub=model.kibria_mix(T,COH_conc,delta_l,w,P_s,P)
df_bub_act=model.kibria_mix_bub_act(T,COH_conc,delta_l,w,P_s,P)
#df_bub=model.kibria_mix_bub_act(T,COH_conc,delta_l,w,P_s,P)

df_bub_tot,ya_bub,yca_bub,y_ohm_bub=model.kibria_mix_bub_act_ohm(T,COH_conc,delta_l,w,P_s,P)

df_comsol=pd.read_csv("COMSOL.csv",index_col=0)
df_exp=pd.read_csv("Experiment.csv",index_col=0)

plt.figure(figsize=(18,9))

plt.plot([i/10 for i in df_comsol.index],df_comsol.iloc[:,0].tolist(),label="COMSOL model (thesis)",linewidth=4,c="r")
plt.plot([i/10 for i in df_exp.index],df_exp.iloc[:,0].tolist(),"o",c='k',markeredgewidth=3, markersize=12,fillstyle="none", label="Experimental data",linewidth=3)  
plt.plot(df_no_bub.index,df_no_bub,label="Kibria Mix Tafel no bubbles activation")
plt.plot(df_bub_act.index,df_bub_act,label="Kibria Mix Tafel with bubbles activation")
plt.plot(df_bub_tot.index,df_bub_tot,label="Kibria Mix Tafel with bubbles activation and ohmic")


plt.legend()
o="$P_{O2}$"
plt.title(f"Set-up: T={T} [K],COH={COH_conc} $[mol/dm^3]$,"+o+f"={P}[atm], w={100*w}[%],$P_s$={P_s} [bar]")
plt.tight_layout()
plt.grid()
#%%



#%%prova

x=np.linspace(0,30000,100)
a=[]
b=[]
T=353
for i in x:
    a.append(model.theta_eps(T,i)[0])
    b.append(model.theta_eps(T,i)[1])
plt.figure(897,figsize=(15,7))    
plt.plot(x,a,label="$\Theta$")
plt.plot(x,b,label="$e$")
plt.legend()

#%% liquid
# a=[]
# for i in np.linspace(70,2000,100):
    
    
#     #a.append(model.sigma_liquid(COH_conc,T)*delta_l*i/1000)
#     a.append(model.eta_liquid(COH_conc,T,i,delta_l))

# plt.figure(5,figsize=(15,8))    
# plt.plot(np.linspace(70,2000,100),a)
# plt.grid()
# plt.xlabel("J $[mA/cm^2]$")
# plt.ylabel("$\eta_l [V]$")    


#%%
# model.summary_kibria_high(T,w,COH_conc,P,P_s)

#%%
# rev=model.v_rev(T,w,P_s)



#%%
# model.eta_sep(COH_conc,T,w,P_s)

#%%
kib_mixed,miles_mix=model.multi_T(w,COH_conc,P,P_s)

plt.figure(figsize=(18,9))
ax=plt.gca()
plt.plot([i/10 for i in df_comsol.index],df_comsol.iloc[:,0].tolist(),label="COMSOL model (thesis)",linewidth=4,c="#ff7f0e")
plt.plot([i/10 for i in df_exp.index],df_exp.iloc[:,0].tolist(),"o",c='k',markeredgewidth=3, markersize=12,fillstyle="none", label="Experimental data",linewidth=3)  
kib_mixed.plot(ax=ax)
ax.set_xlabel("J $[mA/cm^2]$")
ax.set_ylabel("$V_{cell}$ $[V]$")
plt.legend()
plt.grid()
plt.title("Multi temperature")



#%%

kib_100_TA_mix=kib_mixed.loc[100]
Miles_100_TA=miles_mix.loc[100]


dfe=pd.read_csv("Temperature variation overvoltage.csv",index_col=0)

dfe=dfe[dfe.index<=102]


temp=[i-273 for i in np.arange(303,373+5,5)]
re=[]
for i in np.arange(303,373+5,5):
    re.append(model.v_rev(i,w,P_s))
kib_over_mix=kib_100_TA_mix-re
miles_over_mix=Miles_100_TA-re


plt.figure(figsize=(15,8))

plt.plot(temp,[i*1000 for i in kib_over_mix.tolist()],marker="s",label="Kibria Tafel mix")
plt.plot(temp,[i*1000 for i in miles_over_mix.tolist()],marker="s",label="Miles Tafel mix")


plt.plot(dfe.index,dfe.iloc[:,0],label="Experiments data from Henrik",marker="o",linewidth=3,alpha=0.5)

plt.xlabel("Temp [°C]")
plt.ylabel("Overvoltage [mV]")
plt.title("Overvoltage Temp variation @ 0.1 $[A/cm^2]$ with no bubbles")
plt.legend()
plt.grid()



  

    
    
    

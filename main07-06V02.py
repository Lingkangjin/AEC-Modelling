# -*- coding: utf-8 -*-
"""
Created on Thu May 26 23:50:41 2022

@author: utente
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from func_08_06_class import *



#%% using only tafel and semplify  everything
"""
Documentation for this model:
    Intial conditions to be set:
        T: Setting temperature [K]
        P_s: setting system operating pressure [bar]
        w: weight ratio of KOH in electrolyte [0-1]
        delta_l: reaction distance [cm]
    Input:        
        x: a list of the current densities [mA/cm^2]
    output:
        df: a dataframe containing polarisation curve with different conditions
    
    Experiments:
        Kibria experiments performed with: 30 %wt KOH 
        Miles experiments performed with:50 wt% KOH
        
    List of possible outputs:
        -model.kibria(x): Kibria set of i0 and alpha (experiments at 30% wt KOH ) without bubble effects
        -model.kibria_bub_act_ohm(x): Kibria conditions with high overvoltage region including bubble effects
        -model.Miles(x):Miles conditions with high overvoltage region without bubbles
        -model.Miles_bub_act(x): Miles conditions with high overvoltage  work regfion including only activation bubble effects
        -model.Miles_bub_act_ohm(x): Miles conditions with high overvoltage work region including both bubble effects
        -model.kibria_mix_bub_act_ohm(x): Kibria with both sets of i0 and alphas including bubble effects
        -model.kibria_mix_bub_act(x): Kibria with both sets of i0 and alphas including only activation bubble effects
       
    As starting point the normal kibria condition(model.kibria(x)) without the bubble effects should be a good approximation    
"""




#%%folders names and loading experiments data

param="Fitting Parameters" #Parameters fodler name
exp="Experiments"

df_comsol=pd.read_csv(exp+"\\COMSOL.csv",index_col=0)
df_exp=pd.read_csv(exp+"\\Experiment.csv",index_col=0)
dfe=pd.read_csv(exp+"\\Temperature variation overvoltage.csv",index_col=0)

#%% plot configurations
plt.rcParams.update({'font.size': 18, 'font.family': "Times New Roman"})


#%% setting conditions
T=353
# COH_conc=4.212 #mol/dm^3
# P=20 #oxygen pressure 20 atm
P_s=7 #bar system pressure bar
w=0.21
delta_l=0.0132+0.0350 #from thesis (separator and electrode thichness) expressed in [cm]


#%% loading model class
model= AEC_model(T,P_s,w,delta_l)
model.wt_to_cOH()
#%% initialise the current density list


x=np.arange(10,2000,10)

#%% get some polarisation curves



df=model.kibria(x)
df2=model.kibria_bub_act_ohm(x)
df3=model.Miles(x)

#%% comparison plot

plt.figure(figsize=(18,9))

plt.plot([i/10 for i in df_comsol.index],df_comsol.iloc[:,0].tolist(),label="COMSOL model (thesis)",linewidth=4,c="#ff7f0e")
plt.plot([i/10 for i in df_exp.index],df_exp.iloc[:,0].tolist(),"o",c='k',markeredgewidth=3, markersize=12,fillstyle="none", label="Experimental data",linewidth=3)  
plt.plot(df.index,df,label=" Kibria Tafel")
plt.plot(df2.index,df2,label=" Kibria Tafel with Bubbles")
plt.plot(df3.index,df3,label=" Miles Tafel without Bubbles")

plt.legend()
plt.grid()
plt.xlabel("J $[mA/cm^2]$")
plt.ylabel("$V_{cell}$ [V]")
plt.title(f"T={T}[K],w={w*100}%KOH, P={P_s} [bar] delta={round(delta_l,3)} [cm] ")


#%% Teemperature evalutation at 100 mA/cm^2


dfe=dfe[dfe.index<=102]


temp=[i-273 for i in np.arange(303,373+5,5)]
re=[]
   

kib=pd.DataFrame(index=x)
kib_bub=pd.DataFrame(index=x)
miles=pd.DataFrame(index=x)

for j in np.arange(303,373+5,5):
    model= AEC_model(j,P_s,w,delta_l)
    model.wt_to_cOH()
    re.append(model.v_rev())
    kib["temp "+str(j)]=model.kibria(x)
    kib_bub["temp "+str(j)]=model.kibria_bub_act_ohm(x)
    miles["temp "+str(j)]=model.Miles(x)
    
    
    
    
kib_over=kib.loc[100]-re
kib_bub_over=kib_bub.loc[100]-re
miles_over=miles.loc[100]-re


plt.figure(figsize=(15,8))

plt.plot(temp,[i*1000 for i in kib_over.tolist()],marker="s",label="Kibria")
plt.plot(temp,[i*1000 for i in kib_bub_over.tolist()],marker="s",label="Kibria with bubbles")

plt.plot(temp,[i*1000 for i in miles_over.tolist()],marker="s",label="Miles with bubbles")


plt.plot(dfe.index,dfe.iloc[:,0],label="Experiments data from Henrik",marker="o",linewidth=3,alpha=0.5)

plt.xlabel("Temp [°C]")
plt.ylabel("Overvoltage [mV]")
plt.title("Overvoltage Temp variation @ 0.1 $[A/cm^2]$ with no bubbles")
plt.legend()
plt.grid()
  

    
    
    

# -*- coding: utf-8 -*-
"""
Created on Thu May 26 23:50:41 2022

@author: Lingkang Jin l.jin@pm.univpm.it/s221045@dtu.dk
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from func_21_06_class_sakas import *
from Thermal_assess import *



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

    update in 09-06: wt correction in kibria conditions high voltage
    update in 21-06: Activation energy correction
"""




#%%folders names and loading experiments data





param="Fitting Parameters" #Parameters fodler name
exp="Experiments"

df_comsol=pd.read_csv(exp+"\\COMSOL.csv",index_col=0)
df_exp=pd.read_csv(exp+"\\Experiment.csv",index_col=0)
dfe=pd.read_csv(exp+"\\Temperature variation overvoltage.csv",index_col=0)
df_sakas_50=pd.read_csv(exp+"\\Sakas 59.6 °C.csv",index_col=0)
df_sakas_60=pd.read_csv(exp+"\\Sakas 61.15 °C.csv",index_col=0)
df_sakas_70=pd.read_csv(exp+"\\Sakas 70°C .csv",index_col=0)



#%% plot configurations
plt.rcParams.update({'font.size': 18, 'font.family': "Times New Roman"})



#%% sakas setting conditions

T=30+273
# COH_conc=4.212 #mol/dm^3
# P=20 #oxygen pressure 20 atm
P_s=16 #bar system pressure bar
w=0.25
delta_l=0.05+0.08 #sakas
  #sakas

#delta_l=0.0132+0.0350 #from thesis (separator and electrode thichness) expressed in [cm]




Th=Thermal_properties(T)
F=96500 #[C/mol]
v_tn=((Th.delta_h())/(2*F))

#%% loading model class
model= AEC_model(T,P_s,w,delta_l)
model.wt_to_cOH()
#%% initialise the current density list


x=np.linspace(2,500,100)




#%% get some polarisation curves
df=model.kibria_mono(x)


model= AEC_model(59.6+273,P_s,w,delta_l)
model.wt_to_cOH()
df=model.kibria_mono(x)

# df=model.kibria(x)

df_bub_50=model.kibria_bub_ohm(x)

model= AEC_model(61.15+273,P_s,w,delta_l)
model.wt_to_cOH()
df2=model.kibria(x)
df_bub_60=model.kibria_bub_ohm(x)

model= AEC_model(70+273,P_s,w,delta_l)
model.wt_to_cOH()
df3=model.kibria(x)
df_bub_70=model.kibria_bub_ohm(x)


    

#%% comparison plot

plt.figure(figsize=(18,9))


plt.scatter([i*1000 for i in df_sakas_50.index],df_sakas_50.iloc[:,0],marker="o", label="Sakas 59.6°C data")
plt.scatter([i*1000 for i in df_sakas_60.index],df_sakas_60.iloc[:,0],marker="s", label="Sakas 61.15°C data")
plt.scatter([i*1000 for i in df_sakas_70.index],df_sakas_70.iloc[:,0],marker="*", label="Sakas 70°C data")



plt.plot(df_bub_50.index,[i for i in df_bub_50],"--",label=" Kibria Tafel with bubbles 59.6 °C")

plt.plot(df_bub_60.index,[i for i in df_bub_60],"-.",label=" Kibria Tafel with bubbles 61.15 °C")

plt.plot(df_bub_70.index,[i for i in df_bub_70],"-",label=" Kibria Tafel with bubbles 70 °C")

plt.plot(x,[v_tn]*len(x),color="r",linewidth=2,label="Thermal neutral voltage [V]")


plt.ylim(1.25,2.25)

plt.legend()
plt.grid()
plt.xlabel("J $[mA/cm^2]$")
plt.ylabel("$V_{cell}$ [V]")
#T={T}[K]
co="$[mol/dm^3]$"
plt.title(f"w={w*100}%KOH, P={P_s} [bar], delta={round(delta_l,3)} [cm], COH={round(model.COH_conc,1)} "+co+" activation energy correction")



    
    
    

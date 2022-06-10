# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 11:02:25 2022

@author: rafnn
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from func_08_06_class import *


#%% setting conditions
T=65+273.15
P_s=25 #bar system pressure bar
w=0.30
delta_l=(0.0132+0.0350) #from thesis (separator and electrode thichness) expressed in [cm]


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

exp="Experiments"
df_exp=pd.read_csv(exp+"\\ursua.csv",index_col=0)
plt.plot([i*1000 for i in df_exp.index],df_exp.iloc[:,0].tolist(),"o",c='k',markeredgewidth=1, markersize=12,fillstyle="none", label="Experimental data",linewidth=3)  
plt.plot(df.index,df,label=" Kibria Tafel")
plt.plot(df2.index,df2,label=" Kibria Tafel with Bubbles")
plt.plot(df3.index,df3,label=" Miles Tafel without Bubbles")

plt.legend()
plt.grid()
plt.xlabel("J $[mA/cm^2]$")
plt.ylabel("$V_{cell}$ [V]")
#plt.title(f"T={T}[K],w={w*100}%KOH, P={P_s} [bar] delta={round(delta_l,3)} [cm] ")


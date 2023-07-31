# -*- coding: utf-8 -*-
"""
Created on Thu May 26 23:50:41 2022

@author: utente
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from func import *
from Thermal_assess import *

from scipy.optimize import curve_fit

from scipy.optimize import fsolve
efrom sklearn.metrics import r2_score

import seaborn as sns
#%%

plt.rcParams.update({'font.family': "Times New Roman"})
plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.dpi'] = 150

# plt.style.use('classic')
# plt.rcParams.update({'font.family': "Times New Roman"})

# plt.rcParams['font.size'] = 30
# plt.rcParams['legend.fontsize'] = 30
# plt.rcParams['xtick.minor.visible'] = False

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


#%% constants
R=8.31 # [J/mol*K] universal gas constant
F=96500 #[C/mol]  Faraday constant

#%%folders names and loading experiments data





param="Fitting Parameters" #Parameters fodler name
exp="Experiments"

df_sakas_50=pd.read_csv(exp+"\\Sakas 59.6 °C.csv",index_col=0)
df_sakas_60=pd.read_csv(exp+"\\Sakas 61.15 °C.csv",index_col=0)
df_sakas_70=pd.read_csv(exp+"\\Sakas 70°C .csv",index_col=0)

#%% plot configurations
# plt.rcParams.update({'font.size': 18, 'font.family': "Times New Roman"})



#%% sakas setting conditions

T=353
P_s=16 #bar system pressure bar
w=0.25
delta_l=0.05 #sakas
#%% loading thermal properties
Th=Thermal_properties(T)
v_tn=((Th.delta_h())/(2*F))

#%% loading model class
model= AEC_model(T,P_s,w)
model.wt_to_cOH()
#%% initialise the current density list


x=np.linspace(2,350,100)



#%%

Sakas_exp=pd.read_csv(exp+"\\SAKAS.csv")
# Sakas_exp["T [°C]"]=round(Sakas_exp["T [K]"]-273,2)
sampling_number=20
colo=['b', 'g', 'r', 'c', 'm', 'y', 'k']

#%% scatter of experiments
plt.figure(figsize=(10,8),constrained_layout=True)
sns.scatterplot(data=Sakas_exp, x="J [mA/cm^2]", y="V [V]",hue="T [K]",palette=colo,s=80)
plt.grid()
plt.xlabel("$J [mA/cm^2]$")
plt.ylabel("$V_{cell} [V]$")
plt.legend( loc='upper left',title="T [K]")


#%%

# =============================================================================
# Random selecting data, with at leat two different temperatures
# =============================================================================
rand=Sakas_exp.sample(n=sampling_number)

while rand.groupby("T [K]").size().min()<=1 or len(rand["T [K]"].unique())<2:
    rand=Sakas_exp.sample(n=sampling_number)

 
# =============================================================================
# Reversible voltage from model computation
# =============================================================================
        

v_rev=model.v_rev() #V
sigma_l=model.sigma_liquid() #[S/cm]





def func_fit1(X, delta_l, alfa, C1,C2):
    j,T = X
    model= AEC_model(T,P_s,w)
    model.wt_to_cOH()
    sigma_l=model.sigma_liquid()    
    v_rev=model.v_rev() #V

    return v_rev+j*(delta_l/(sigma_l*1000))+((R*T)/(alfa*2*F))*(np.log(j)-(C1/T+C2))


T_s=Sakas_exp["T [K]"].unique()




j=rand.iloc[:,0].to_numpy()
T=rand.iloc[:,1].to_numpy()
y=rand.iloc[:,-1].to_numpy()

# =============================================================================
# # intitial guess and boundaries
# =============================================================================
p0=0.40,0.13,-6330,14  #delta_l, alfa, C1,C2
bnds=((0,0.1,-8000,0),(1,0.4,-3000,15))

# =============================================================================
# Finding coefficients
# =============================================================================
popt, pcov=curve_fit(func_fit1, (j,T), y, p0,bounds=bnds)

j_new=np.linspace(50,400,100)

#%%

plt.style.use('classic')

plt.rcParams['legend.fontsize'] = 16
plt.rcParams['xtick.minor.visible'] = False
plt.rcParams.update({'font.size': 18, 'font.family': "Times New Roman"})



plt.figure(figsize=(6,5),constrained_layout=True)
sns.scatterplot(data=Sakas_exp, x="J [mA/cm^2]", y="V [V]",hue="T [K]",palette=colo,legend=False)
sns.scatterplot(data=rand, x="J [mA/cm^2]", y="V [V]",hue="T [K]",marker="s",s=80,palette=colo)

ax=plt.gca()

for l in T_s:
    T_new=[l]*len(j_new)
    
    y_fit=[]
    
    for i in range(len(j_new)):        
        y_fit.append(func_fit1((j_new[i],T_new[i]),*popt))
        
    d=Sakas_exp[Sakas_exp["T [K]"]==l]
    
    y_pred=[]
    for jj in d["J [mA/cm^2]"]:
        y_pred.append(func_fit1((jj,l),*popt))
        
    y_real=d["V [V]"]    
    
    s=r2_score(y_real, y_pred)    


     
    ax.plot(j_new,y_fit,label=f"fit {round(l,2)} [K]")
    # ax.plot(j_new,y_fit,label=f"fitted @ {l} [K],$R_2$={round(s,3)}")


ax.legend()
ax.set_title("Sakas")


ax.grid()


# Put a legend to the right of the current axis
ax.legend(ncol=2,loc=0)
ax.set_xlabel("$J$ $[mA/cm^2]$")
ax.set_ylabel("$V_{cell}$ $[V]$")


print(rand)


#%%

df_univpm=pd.read_csv(exp+"\\UNIVPM_data.csv",index_col=0)




#%%

df=df_univpm.copy()

# cols=['AIM_Stack_1_temp','Temp. cat', 'Cell voltage [V]','J from faraday $[mA/cm^2]$','AIM_Gas_pressure_B2'
#       ,'Pres. cat [bar]','AIM_Gas_pressure_B2']

# for i in df.columns:
#     if i in cols:
#         pass
#     else:
#         df.drop(columns=[i],axis=1, inplace=True)
        
        
df=df[df['Pres. cat [bar]']==4.5]        #4 bar

#%%
df["Temp. cat2"] = pd.cut(
    x=df["AIM_Stack_1_temp"],
    bins=np.arange(-5, 60, 5),
    labels=np.arange(0, 60, 5),
)

df=df[df["Temp. cat"]==50]
#%%

df["T [K]"]=df['Temp. cat']+273.15
#%%
df_test=df[['J from faraday $[mA/cm^2]$',"T [K]",'Cell voltage [V]','Temp. cat']]        

#%%
plt.figure()
sns.scatterplot(data=df_test, x="J from faraday $[mA/cm^2]$", y="Cell voltage [V]",hue="Temp. cat",
                palette=colo,s=80)
#%%
plt.scatter(df_test[df_test['Temp. cat']==50]['J from faraday $[mA/cm^2]$'],df_test[df_test['Temp. cat']==50]['Cell voltage [V]'])

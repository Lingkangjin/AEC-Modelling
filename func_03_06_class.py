# -*- coding: utf-8 -*-
"""
Created on Tue May 24 15:39:53 2022

@author: utente
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from sklearn.metrics import r2_score
import pickle



#%%
#solving io with arrenhius plot for tafel

# Kibria with mixed io

#Bubble effects 





#%%
df_c=pd.read_pickle("df_c.pkl")
DF_ARR=pd.read_pickle("DF_ARR.pkl")

# with open("i0.txt", "rb") as myFile:
#     io_par = pickle.load(myFile)



with open("tafel.txt", "rb") as myFile:
    tafel_par = pickle.load(myFile)
    
    

with open("alpha.txt", "rb") as myFile:
    alpha_par = pickle.load(myFile)    

#%%

# T=353
# COH_conc=4.212 #mol/dm^3
# P=20 #oxygen pressure 20 atm
# P_s=7 #bar system pressure bar
# w=0.35
#%%


class AEC_model(object):

# T=353
# COH_conc=4.212 #mol/dm^3
# P=20 #oxygen pressure 20 atm
# P_s=7 #bar system pressure bar
# w=0.35
# delta_l=0.0132+0.0350 #from thesis (separator and electrode thichness)    

    def __init__(self, T,COH_conc,P,P_s,w,delta_l):
        self.T=T
        self.COH_conc=COH_conc
        self.P=P
        self.P_s=P_s
        self.w=w
        self.delta_l=delta_l
    #%% bubble surface coverage 

    def theta_eps(self,T,J):
        J_lim=30000   # Piontelli results [mA/cm^2]
        T_0=273
        #theta=(-97.25+182*(T/T_0)-84*(T/T_0)**2)*(J/J_lim)**(0.3) #To be revised temperature dependence T_0
        theta=0.045*(J)**(0.3)
        #theta=0.023*(J/J_lim)**(0.3)
        eps=theta*(2/3)
        return theta,eps
    #%% current density limit
    
    def anode_i_limit(self,T,df_c):
        key="OER"    
    
        return np.polyval(df_c[key],T)
    
    
    def cathode_i_limit(self,T,df_c):
        key="HER"
        df_c=pd.read_pickle("df_c.pkl")
    
        return np.polyval(df_c[key],T)
    #%% arrenhius exchange current density
    def io_func(self,key,T):
        T_ref=DF_ARR.loc[key][0]
        io_ref=DF_ARR.loc[key][1]
        R=0.00831 #[kJ/mol K]
        E=DF_ARR.loc[key][2]
        
        return io_ref*np.exp((E/R)*((1/T_ref)-(1/T)))
   #%% reversible voltage calculation     
    def reve_V(self,T_in):
            rev=1.5184-1.5421*10**(-3)*T_in+9.523*10**(-5)*T_in*np.log(T_in)+9.84*10**(-8)*T_in**2
            return rev
        
    def v_rev(self,T_in,w,P_s):
            weight_conc=w #[kgKoH/kgl]
            KOH_mol_w=56.1        #[g/mol]
            rho=self.density(self.A(T_in-273),w)   #A function is a polynomial using Temperature in Celsius
            m=weight_conc*rho/KOH_mol_w
            a=-0.0151*m-1.6788*10**(-3)*m**2+2.2588*10**(-5)*m**3
            b=1-1.2062*10**(-3)*m+5.6024*10**(-4)*m**2-7.8228*10**(-6)*m**(3)
            Pv_h20=np.exp(81.6179-7699.68/T_in-10.9*np.log(T_in)+9.5891*10**(-3)*T_in)
            Pv_KOH=np.exp(2.302*a+b*np.log(Pv_h20))
            alfa_H20=np.exp(-0.05192*m+0.003302*m**2+(3.177*m-2.131*m**2)/T_in)
            R=8.31 # [J/mol*K]
            F=96500 #[C/mol]
            # P=P_s #bar
            V_re=self.reve_V(T_in)+(R*T_in)/(2*F)*np.log((P_s-Pv_KOH)**(1.5)/alfa_H20)
            return V_re  
        
    #%% plot
    
    plt.rcParams.update({'font.size': 18, 'font.family': "Times New Roman"})
    
    # def plot_rev(self,i,w,P_s):
    #     T_range=[i for i in range(298,350)]
        
    #     V_rev_range=[]
    #     for i in T_range:    
    #         V_rev_range.append(self.v_rev(i,w,P_s))
        
        
    #     plt.figure(figsize=(8,5))
    #     plt.plot(T_range,V_rev_range,marker='.')
    #     plt.xlabel("Temp [K]")
    #     plt.ylabel("$V_{rev}$")
    #     plt.grid()    
        
    #%% anodic overpotential
    
    
    
    def anodic_over_Kibria_low(self,T_in,J):#J expressed in [mA/cm^2]   
    
        key= 'Kibria anode-low'
        i_o=self.io_func(key,T_in)
        A=np.polyval(tafel_par[key],T_in) 
        eta_an=A*np.log10(J/i_o)
        return eta_an,i_o
    
    def anodic_over_Kibria_high(self,T_in,J):#J expressed in [mA/cm^2] 
        key= 'Kibria anode-high'
        i_o=self.io_func(key,T_in)
        A=np.polyval(tafel_par[key],T_in)         
        eta_an=A*np.log10(J/i_o)
        return eta_an,i_o
    
    def anodic_over_miles_low(self,T_in,J):#J expressed in [mA/cm^2]
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]
        key='Miles anode-low'
        i_o=self.io_func(key,T_in) #Arrenhius
        # i_o=np.polyval(io_par[key],T_in) #polynomial
        alpha=np.polyval(alpha_par[key],T_in)    
        
        eta_an=2.3*R*T_in/(alpha*F)*np.log10(J/i_o)
        return eta_an,i_o    
    
    
    
    
    def anodic_over_miles_high(self,T_in,J):#J expressed in [mA/cm^2]
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]
        key='Miles anode-high'
        i_o=self.io_func(key,T_in)
        
        #i_o=np.polyval(io_par[key],T_in) 
        alpha=np.polyval(alpha_par[key],T_in)
        eta_an=2.3*R*T_in/(alpha*F)*np.log10(J/i_o)
        return eta_an,i_o    
    #%% anodic bubbles overpotential
    
    
    
    def anodic_over_Bubbles_low(self,T_in,J):#J expressed in [mA/cm^2]   
    
        key= 'Kibria anode-low'
        i_o=self.io_func(key,T_in)
        A=np.polyval(tafel_par[key],T_in) 
        Jeff=J/(1-self.theta_eps(T_in,J)[0])
        eta_an_bu=A*np.log10(Jeff/i_o)
        return eta_an_bu
    
    def anodic_over_Bubbles_high(self,T_in,J):#J expressed in [mA/cm^2] 
        key= 'Kibria anode-high'
        A=np.polyval(tafel_par[key],T_in)     
        i_o=self.io_func(key,T_in)
        Jeff=J/(1-self.theta_eps(T_in,J)[0])
        eta_an_bu=A*np.log10(Jeff/i_o)
        return eta_an_bu
    
    
    def anodic_over_Bubbles_low_miles(self,T_in,J):#J expressed in [mA/cm^2]   
    
        key= 'Miles anode-low'
        i_o=self.io_func(key,T_in)
        A=np.polyval(tafel_par[key],T_in) 
        Jeff=J/(1-self.theta_eps(T_in,J)[0])
        eta_an_bu=A*np.log10(Jeff/i_o)
        return eta_an_bu
    
    def anodic_over_Bubbles_high_miles(self,T_in,J):#J expressed in [mA/cm^2] 
        key= 'Miles anode-high'
        A=np.polyval(tafel_par[key],T_in)     
        i_o=self.io_func(key,T_in)
        Jeff=J/(1-self.theta_eps(T_in,J)[0])
        eta_an_bu=A*np.log10(Jeff/i_o)
        return eta_an_bu
   
    #%% anodic overpotential plot
    # def an_plot(self,T):
    #     x=np.linspace(70,2000,100)
    #     y=[]
    #     y1=[]
    #     y2=[]
    #     y3=[]
    #     for i in x:       
    #         y.append(self.anodic_over_Kibria_low(T,i)[0])
    #         y1.append(self.anodic_over_Kibria_high(T,i)[0])
    #         y2.append(self.anodic_over_miles_low(T,i)[0])
    #         y3.append(self.anodic_over_miles_high(T,i)[0])
         
        
    #     fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(16,7))
        
        
    #     ax1.plot(x,y,marker='.',label="low $\eta$")
    #     ax1.plot(x,y1,marker=".",label="high $\eta$")
    #     ax1.set_xlabel("J $[mA/cm^2]$")
    #     ax1.set_ylabel("$\eta_{an}$")
    #     ax1.legend()
    #     ax1.set_title("Kibria, wt(KOH)=30%")
    #     ax1.grid()
        
    #     ax2.plot(x,y2,marker='.',label="low $\eta$")
    #     ax2.plot(x,y3,marker='.',label="high $\eta$")
    #     ax2.set_xlabel("J $[mA/cm^2]$")
    #     ax2.set_ylabel("$\eta_{an} [V]$")
    #     ax2.grid()
    #     ax2.set_title("Miles, wt(KOH)=50%")
    #     ax2.legend()
        
    #     plt.suptitle("Anodic overpotential")
    #%% cathodic over potnetial
    
    def cathodic_over_Kibria_low(self,T_in,J):#J expressed in [mA/cm^2]   
        key= 'Kibria cathode-low'
        i_o=self.io_func(key,T_in)
        # i_o=np.polyval(io_par[key],T_in)    
        A=np.polyval(tafel_par[key],T_in)    
        eta_an=A*np.log10(J/i_o)
        return eta_an,i_o
    
    def cathodic_over_Kibria_high(self,T_in,J):#J expressed in [mA/cm^2]   
        key= 'Kibria cathode-high'
        i_o=self.io_func(key,T_in)
        #i_o=np.polyval(io_par[key],T_in)    
        A=np.polyval(tafel_par[key],T_in)    
        eta_an=A*np.log10(J/i_o)
        return eta_an,i_o
    
    def cathodic_over_miles_low(self,T_in,J):#J expressed in [mA/cm^2]
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]
        
        key='Miles cathode-low'
        i_o=self.io_func(key,T_in)
        # i_o=np.polyval(io_par[key],T_in) 
        alpha=np.polyval(alpha_par[key],T_in)    
        
        eta_ca=2.3*R*T_in/(alpha*F)*np.log10(J/i_o)
        return eta_ca,i_o
    
    
    def cathodic_over_miles_high(self,T_in,J):#J expressed in [mA/cm^2]
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]
        
        key='Miles cathode-high'
        i_o=self.io_func(key,T_in)
        # i_o=np.polyval(io_par[key],T_in) 
        alpha=np.polyval(alpha_par[key],T_in)    
        
       
        eta_ca=2.3*R*T_in/(alpha*F)*np.log10(J/i_o)
        return eta_ca,i_o
    #%% cathodic overpotential due to bubbles
    
    def cathodic_over_Bubbles_low(self,T_in,J):#J expressed in [mA/cm^2]   
        key= 'Kibria cathode-low'
        A=np.polyval(tafel_par[key],T_in)
        i_o=self.io_func(key,T_in)
        Jeff=J/(1-self.theta_eps(T_in,J)[0])
        eta_an_bub=A*np.log10(Jeff/i_o)
        return eta_an_bub
    
    def cathodic_over_Bubbles_high(self,T_in,J):#J expressed in [mA/cm^2]   
        key= 'Kibria cathode-high'
        i_o=self.io_func(key,T_in)
        A=np.polyval(tafel_par[key],T_in)    
        Jeff=J/(1-self.theta_eps(T_in,J)[0])
        eta_an_bub=A*np.log10(Jeff/i_o)
        return eta_an_bub
    
    def cathodic_over_Bubbles_low_miles(self,T_in,J):#J expressed in [mA/cm^2]   
        key= 'Miles cathode-low'
        A=np.polyval(tafel_par[key],T_in)
        i_o=self.io_func(key,T_in)
        Jeff=J/(1-self.theta_eps(T_in,J)[0])
        eta_an_bub=A*np.log10(Jeff/i_o)
        return eta_an_bub
    
    def cathodic_over_Bubbles_high_miles(self,T_in,J):#J expressed in [mA/cm^2]   
        key= 'Miles cathode-high'
        i_o=self.io_func(key,T_in)
        A=np.polyval(tafel_par[key],T_in)    
        Jeff=J/(1-self.theta_eps(T_in,J)[0])
        eta_an_bub=A*np.log10(Jeff/i_o)
        return eta_an_bub
    
    
    #%% cathodic over potential  plot
    
    # def ca_plot(self,T):
    #     x=np.arange(70,500,10)
    #     y=[]
    #     y1=[]
    #     y2=[]
    #     y3=[]
    #     for i in x:        
    #         y.append(self.cathodic_over_Kibria_low(T,i)[0])
    #         y1.append(self.cathodic_over_Kibria_high(T,i)[0])
    #         y2.append(self.cathodic_over_miles_low(T,i)[0])
    #         y3.append(self.cathodic_over_miles_high(T,i)[0])
         
        
    #     fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(16,7))
        
        
    #     ax1.plot(x,y,marker='.',label="low $\eta$")
    #     ax1.plot(x,y1,marker=".",label="high $\eta$")
    #     ax1.set_xlabel("J $[mA/cm^2]$")
    #     ax1.set_ylabel("$\eta_{ca}$")
    #     ax1.legend()
    #     ax1.set_title("Kibria, wt(KOH)=30%")
    #     ax1.grid()
        
    #     ax2.plot(x,y2,marker='.',label="low $\eta$")
    #     ax2.plot(x,y3,marker='.',label="high $\eta$")
    #     ax2.set_xlabel("J $[mA/cm^2]$")
    #     ax2.set_ylabel("$\eta_{ca} [V]$")
    #     ax2.legend()
    #     ax2.grid()
    #     ax2.set_title("Miles, wt(KOH)=50%")
        
    #     plt.suptitle("Cathodic overpotential")
    
    #%% Nickel electrode overpotential
    def eta_ni(self,T,delta_l,J):
        return J*delta_l/(1.477*(T-273)**2+791.11*(T-273)+160285)
    #%% Bubbles ohmic resistance
    
    def eta_ohm_bub(self,con,T,J,delta_l): #con in [mol/dm^3] 
        sigma_l_f=self.sigma_liquid(con,T)
        sigma_b_f=sigma_l_f/((1/(1-self.theta_eps(T,J)[1])**(1.5))-1) #Bruggman relationship [S/cm]
        return J*delta_l/(sigma_b_f*1000)
    
    
    
    
    #%% electrolyte ion mass transport overpotential
    
    
    def sigma_liquid(self,con,T_kel): #con=concentration of OH in [mol/dm^3]
        sig_l=-2.041*con-0.0028*con**2+0.005332 *con*T_kel+207.2*con/T_kel+0.001043*con**3-0.0000003*con**2*T_kel**2
        return sig_l #[S/cm] or [A/V*cm]
    
    
   
    
  
        
    def eta_liquid(self,COH_conc,T,J,delta_l): #J in mA/cm^2
        return(delta_l*J/(self.sigma_liquid(COH_conc,T)*1000))
    
    
    #%%  electrolyte overpotential plot 
    def l_plot(self,COH_conc,T,delta_l):
        x=np.linspace(0,2000,100) #mA/cm^2 to be converted in A/cm^2
        y_l=[]
        for i in x:
            y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
        
        plt.figure(5,figsize=(10,7))
        plt.plot(x,y_l,marker=".")
        plt.grid()
        plt.xlabel("$J [mA/cm^2]$")
        plt.ylabel("$\eta_l [V]$")
    
    #%% ION mass tranport through separator  over voltage eta d
    
    
    # def eta_d(self,T,J): #temperature dependence
    #     R=0.060+80*np.exp(T/50)/10000 #0.17 from thesis [ohm * cm^2]
        
    #     eta=R*J/1000 #mA/cm^2 into A/cm^2
    #     self.R_eta_d=R
    #     return eta
    
    def eta_d_Vermeiren(self,J): #VERmeiren model
        R=0.11 #0.17 from thesis [ohm * cm^2]
        eta=R*J/1000 #mA/cm^2 into A/cm^2
        self.R_eta_d=R
        return eta
    
    
    #%% #Overvoltage  evalutation
    
    # def d_plot(self,):
    #     x=np.linspace(70,2000,100) #mA/cm^2 to be converted in A/cm^2
    #     y_d=[]
    #     for i in x:
    #         y_d.append(self.eta_d_Vermeiren(i))
        
    #     plt.figure(5,figsize=(10,7))
    #     plt.plot(x,y_d,marker=".")
    #     plt.grid()
    #     plt.xlabel("$J [mA/cm^2]$")
    #     plt.ylabel("$\eta_d [V]$")
    
    #%% miles mix
    
    def Miles_mix_bub_act(self,T,COH_conc,delta_l,w,P_s,P): #only high
        x=np.arange(10,2000,10) #mA/cm^2 to be converted in A/cm^2
        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        tot=[]
        an_i_lim=self.anode_i_limit(T,df_c)
        ca_i_lim=self.cathode_i_limit(T,df_c)        
        
        for i in x:
            # if i < an_i_lim and i< ca_i_lim:
            #     ya.append(self.anodic_over_miles_low(T,i)[0])
            #     yca.append(self.cathodic_over_miles_low(T,i)[0])
            #     ya_bub.append(self.anodic_over_Bubbles_low_miles(T,i))
            #     yca_bub.append(self.cathodic_over_Bubbles_low_miles(T,i))
                
            #     y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            #     y_d.append(self.eta_d_Vermeiren(i))
            #     y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
            #     tot.append(self.v_rev(T,w,P_s)+
            #                ya[-1]+
            #                yca[-1]+
            #                y_l[-1]+
            #                y_d[-1])
              
            # elif i < an_i_lim and i> ca_i_lim:
            #     ya.append(self.anodic_over_miles_low(T,i)[0])
            #     yca.append(self.cathodic_over_miles_high(T,i)[0])
            #     ya_bub.append(self.anodic_over_Bubbles_low_miles(T,i))
            #     yca_bub.append(self.cathodic_over_Bubbles_high_miles(T,i))
            #     y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            #     y_d.append(self.eta_d_Vermeiren(i))
            #     y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
            #     tot.append(self.v_rev(T,w,P_s)+
            #                ya[-1]+
            #                yca[-1]+
            #                y_l[-1]+
            #                y_d[-1])
            # elif i > an_i_lim and i< ca_i_lim:
            #     ya.append(self.anodic_over_miles_high(T,i)[0])
            #     yca.append(self.cathodic_over_miles_low(T,i)[0])
            #     ya_bub.append(self.anodic_over_Bubbles_high_miles(T,i))
            #     yca_bub.append(self.cathodic_over_Bubbles_low_miles(T,i))
            #     y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
            #     y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            #     y_d.append(self.eta_d_Vermeiren(i))
                
            #     tot.append(self.v_rev(T,w,P_s)+
            #                ya[-1]+
            #                yca[-1]+
            #                y_l[-1]+
            #                y_d[-1])
                
            # elif i > an_i_lim and i> ca_i_lim:
                # ya.append(self.anodic_over_miles_high(T,i)[0])
                # yca.append(self.cathodic_over_miles_high(T,i)[0])
                # ya_bub.append(self.anodic_over_Bubbles_high_miles(T,i))
                # yca_bub.append(self.cathodic_over_Bubbles_high_miles(T,i))
                # y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
                # y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                # y_d.append(self.eta_d_Vermeiren(i))
                
                # tot.append(self.v_rev(T,w,P_s)+
                #            ya[-1]+
                #            yca[-1]+
                #            y_l[-1]+
                #            y_d[-1])
                # else:
                #     pass
            ya.append(self.anodic_over_miles_high(T,i)[0])
            yca.append(self.cathodic_over_miles_high(T,i)[0])
            ya_bub.append(self.anodic_over_Bubbles_high_miles(T,i))
            yca_bub.append(self.cathodic_over_Bubbles_high_miles(T,i))
            y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
            
            y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            y_d.append(self.eta_d_Vermeiren(i))
            
            tot.append(self.v_rev(T,w,P_s)+
                       ya[-1]+
                       yca[-1]+
                       y_l[-1]+
                       y_d[-1])
            
           
            
                
                
                
                
          
        
        df=pd.DataFrame(index=x)
        
       
        
        df["Tafel"]=tot
        return df["Tafel"]
    
    #%% Miles ohmic
    
    def Miles_mix_bub_ohm(self,T,COH_conc,delta_l,w,P_s,P):
        x=np.arange(10,2000,10) #mA/cm^2 to be converted in A/cm^2
        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        tot=[]
        an_i_lim=self.anode_i_limit(T,df_c)
        ca_i_lim=self.cathode_i_limit(T,df_c)        
        
        for i in x:
            # if i < an_i_lim and i< ca_i_lim:
            #     ya.append(self.anodic_over_miles_low(T,i)[0])
            #     yca.append(self.cathodic_over_miles_low(T,i)[0])
            #     ya_bub.append(self.anodic_over_Bubbles_low_miles(T,i))
            #     yca_bub.append(self.cathodic_over_Bubbles_low_miles(T,i))
                
            #     y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            #     y_d.append(self.eta_d_Vermeiren(i))
            #     y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
            #     tot.append(self.v_rev(T,w,P_s)+
            #                ya[-1]+
            #                yca[-1]+
            #                y_l[-1]+
            #                y_d[-1]+
            #                y_ohm_bub[-1])
              
            # elif i < an_i_lim and i> ca_i_lim:
            #     ya.append(self.anodic_over_miles_low(T,i)[0])
            #     yca.append(self.cathodic_over_miles_high(T,i)[0])
            #     ya_bub.append(self.anodic_over_Bubbles_low_miles(T,i))
            #     yca_bub.append(self.cathodic_over_Bubbles_high_miles(T,i))
            #     y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            #     y_d.append(self.eta_d_Vermeiren(i))
            #     y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
            #     tot.append(self.v_rev(T,w,P_s)+
            #                ya[-1]+
            #                yca[-1]+
            #                y_l[-1]+
            #                y_d[-1]+
            #                y_ohm_bub[-1])
            # elif i > an_i_lim and i< ca_i_lim:
            #     ya.append(self.anodic_over_miles_high(T,i)[0])
            #     yca.append(self.cathodic_over_miles_low(T,i)[0])
            #     ya_bub.append(self.anodic_over_Bubbles_high_miles(T,i))
            #     yca_bub.append(self.cathodic_over_Bubbles_low_miles(T,i))
            #     y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
            #     y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            #     y_d.append(self.eta_d_Vermeiren(i))
                
            #     tot.append(self.v_rev(T,w,P_s)+
            #                ya[-1]+
            #                yca[-1]+
            #                y_l[-1]+
            #                y_d[-1]+
            #                y_ohm_bub[-1])
                
            # elif i > an_i_lim and i> ca_i_lim:
            #     ya.append(self.anodic_over_miles_high(T,i)[0])
            #     yca.append(self.cathodic_over_miles_high(T,i)[0])
            #     ya_bub.append(self.anodic_over_Bubbles_high_miles(T,i))
            #     yca_bub.append(self.cathodic_over_Bubbles_high_miles(T,i))
            #     y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
            #     y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            #     y_d.append(self.eta_d_Vermeiren(i))
                
            #     tot.append(self.v_rev(T,w,P_s)+
            #                ya[-1]+
            #                yca[-1]+
            #                y_l[-1]+
            #                y_d[-1]+
            #                y_ohm_bub[-1])
                
           
            # else:
            #     pass
            ya.append(self.anodic_over_miles_high(T,i)[0])
            yca.append(self.cathodic_over_miles_high(T,i)[0])
            ya_bub.append(self.anodic_over_Bubbles_high_miles(T,i))
            yca_bub.append(self.cathodic_over_Bubbles_high_miles(T,i))
            y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
            
            y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            y_d.append(self.eta_d_Vermeiren(i))
            
            tot.append(self.v_rev(T,w,P_s)+
                       ya[-1]+
                       yca[-1]+
                       y_l[-1]+
                       y_d[-1]+
                       y_ohm_bub[-1])
                
                
                
                
          
        
        df=pd.DataFrame(index=x)
        
       
        
        df["Tafel"]=tot
        return df["Tafel"]
    
    
    #%% Kibria Mix bubbles activation
    def kibria_mix_bub_act(self,T,COH_conc,delta_l,w,P_s,P):
        x=np.arange(10,2000,10) #mA/cm^2 to be converted in A/cm^2
        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        tot=[]
        an_i_lim=self.anode_i_limit(T,df_c)
        ca_i_lim=self.cathode_i_limit(T,df_c)        
        
        for i in x:
            if i < an_i_lim and i< ca_i_lim:
                ya.append(self.anodic_over_Kibria_low(T,i)[0])
                yca.append(self.cathodic_over_Kibria_low(T,i)[0])
                ya_bub.append(self.anodic_over_Bubbles_low(T,i))
                yca_bub.append(self.cathodic_over_Bubbles_low(T,i))
                
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
                tot.append(self.v_rev(T,w,P_s)+
                           #self.anodic_over_Kibria_low(T,i)[0]+
                           self.anodic_over_Bubbles_low(T,i)+
                           #self.cathodic_over_Kibria_low(T,i)[0]+
                           self.cathodic_over_Bubbles_low(T,i)+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                           #self.eta_ohm_bub(COH_conc,T,i,delta_l)+
                           self.eta_ni(T,delta_l,i))
                
                # tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_low(T,i)[0]+
                #            # self.anodic_over_Bubbles_low(T,i)+
                #            self.cathodic_over_Kibria_low(T,i)[0]+
                #            # self.cathodic_over_Bubbles_low(T,i)+
                #            self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                #            self.eta_ohm_bub(COH_conc,T,i,delta_l)+
                #            self.eta_ni(T,delta_l,i))
            elif i < an_i_lim and i> ca_i_lim:
                ya.append(self.anodic_over_Kibria_low(T,i)[0])
                yca.append(self.cathodic_over_Kibria_high(T,i)[0])
                ya_bub.append(self.anodic_over_Bubbles_low(T,i))
                yca_bub.append(self.cathodic_over_Bubbles_high(T,i))
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
                tot.append(self.v_rev(T,w,P_s)+
                           #self.anodic_over_Kibria_low(T,i)[0]+
                           self.anodic_over_Bubbles_low(T,i)+
                           #self.cathodic_over_Kibria_high(T,i)[0]+
                           self.cathodic_over_Bubbles_high(T,i)+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                           #self.eta_ohm_bub(COH_conc,T,i,delta_l)
                           +self.eta_ni(T,delta_l,i))
            elif i > an_i_lim and i< ca_i_lim:
                ya.append(self.anodic_over_Kibria_high(T,i)[0])
                yca.append(self.cathodic_over_Kibria_low(T,i)[0])
                ya_bub.append(self.anodic_over_Bubbles_high(T,i))
                yca_bub.append(self.cathodic_over_Bubbles_low(T,i))
                y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                tot.append(self.v_rev(T,w,P_s)+
                           #self.anodic_over_Kibria_high(T,i)[0]+
                           self.anodic_over_Bubbles_high(T,i)+
                           #self.cathodic_over_Kibria_low(T,i)[0]+
                           self.cathodic_over_Bubbles_low(T,i)+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                           #self.eta_ohm_bub(COH_conc,T,i,delta_l)+
                           self.eta_ni(T,delta_l,i))
            elif i > an_i_lim and i> ca_i_lim:
                ya.append(self.anodic_over_Kibria_high(T,i)[0])
                yca.append(self.cathodic_over_Kibria_high(T,i)[0])
                ya_bub.append(self.anodic_over_Bubbles_high(T,i))
                yca_bub.append(self.cathodic_over_Bubbles_high(T,i))
                y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                tot.append(self.v_rev(T,w,P_s)+
                           #self.anodic_over_Kibria_high(T,i)[0]+
                           self.anodic_over_Bubbles_high(T,i)+
                           #self.cathodic_over_Kibria_high(T,i)[0]+
                           self.cathodic_over_Bubbles_high(T,i)+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                           #self.eta_ohm_bub(COH_conc,T,i,delta_l)+
                           self.eta_ni(T,delta_l,i))
            else:
                pass
                
                
                
                
          
        
        df=pd.DataFrame(index=x)
        
       
        
        df["Tafel"]=tot
        return df["Tafel"]
    
    #%% bubbles activation+ ohmic
    
    
    def kibria_mix_bub_act_ohm(self,T,COH_conc,delta_l,w,P_s,P):
        x=np.arange(10,2000,10) #mA/cm^2 to be converted in A/cm^2
        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        tot=[]
        an_i_lim=self.anode_i_limit(T,df_c)
        ca_i_lim=self.cathode_i_limit(T,df_c)        
        
        for i in x:
            if i < an_i_lim and i< ca_i_lim:
                ya.append(self.anodic_over_Kibria_low(T,i)[0])
                yca.append(self.cathodic_over_Kibria_low(T,i)[0])
                ya_bub.append(self.anodic_over_Bubbles_low(T,i))
                yca_bub.append(self.cathodic_over_Bubbles_low(T,i))
                
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
                tot.append(self.v_rev(T,w,P_s)+
                           #self.anodic_over_Kibria_low(T,i)[0]+
                           self.anodic_over_Bubbles_low(T,i)+
                           #self.cathodic_over_Kibria_low(T,i)[0]+
                           self.cathodic_over_Bubbles_low(T,i)+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                           self.eta_ohm_bub(COH_conc,T,i,delta_l)+
                           self.eta_ni(T,delta_l,i))
                
                # tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_low(T,i)[0]+
                #            # self.anodic_over_Bubbles_low(T,i)+
                #            self.cathodic_over_Kibria_low(T,i)[0]+
                #            # self.cathodic_over_Bubbles_low(T,i)+
                #            self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                #            self.eta_ohm_bub(COH_conc,T,i,delta_l)+
                #            self.eta_ni(T,delta_l,i))
            elif i < an_i_lim and i> ca_i_lim:
                ya.append(self.anodic_over_Kibria_low(T,i)[0])
                yca.append(self.cathodic_over_Kibria_high(T,i)[0])
                ya_bub.append(self.anodic_over_Bubbles_low(T,i))
                yca_bub.append(self.cathodic_over_Bubbles_high(T,i))
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
                tot.append(self.v_rev(T,w,P_s)+
                           #self.anodic_over_Kibria_low(T,i)[0]+
                           self.anodic_over_Bubbles_low(T,i)+
                           #self.cathodic_over_Kibria_high(T,i)[0]+
                           self.cathodic_over_Bubbles_high(T,i)+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                           self.eta_ohm_bub(COH_conc,T,i,delta_l)
                           +self.eta_ni(T,delta_l,i))
            elif i > an_i_lim and i< ca_i_lim:
                ya.append(self.anodic_over_Kibria_high(T,i)[0])
                yca.append(self.cathodic_over_Kibria_low(T,i)[0])
                ya_bub.append(self.anodic_over_Bubbles_high(T,i))
                yca_bub.append(self.cathodic_over_Bubbles_low(T,i))
                y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                tot.append(self.v_rev(T,w,P_s)+
                           #self.anodic_over_Kibria_high(T,i)[0]+
                           self.anodic_over_Bubbles_high(T,i)+
                           #self.cathodic_over_Kibria_low(T,i)[0]+
                           self.cathodic_over_Bubbles_low(T,i)+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                           self.eta_ohm_bub(COH_conc,T,i,delta_l)+
                           self.eta_ni(T,delta_l,i))
            elif i > an_i_lim and i> ca_i_lim:
                ya.append(self.anodic_over_Kibria_high(T,i)[0])
                yca.append(self.cathodic_over_Kibria_high(T,i)[0])
                ya_bub.append(self.anodic_over_Bubbles_high(T,i))
                yca_bub.append(self.cathodic_over_Bubbles_high(T,i))
                y_ohm_bub.append(self.eta_ohm_bub(COH_conc,T,i,delta_l))
                
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                tot.append(self.v_rev(T,w,P_s)+
                           #self.anodic_over_Kibria_high(T,i)[0]+
                           self.anodic_over_Bubbles_high(T,i)+
                           #self.cathodic_over_Kibria_high(T,i)[0]+
                           self.cathodic_over_Bubbles_high(T,i)+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i)+
                           self.eta_ohm_bub(COH_conc,T,i,delta_l)+
                           self.eta_ni(T,delta_l,i))
            else:
                pass
                
                
                
                
          
        
        df=pd.DataFrame(index=x)
        
       
        
        df["Tafel"]=tot
        
       
        
        
        return df["Tafel"],ya_bub,yca_bub,y_ohm_bub
    
    
    
    #%% bubble effects on activation
    # def kibria_mix_bub_act(self,T,COH_conc,delta_l,w,P_s,P):
    #     x=np.arange(10,2000,10) #mA/cm^2 
    #     ya=[]
    #     yca=[]
    #     ya_bub=[]
    #     yca_bub=[]
        
    #     y_l=[]
    #     y_d=[]
    #     tot=[]
    #     an_i_lim=self.anode_i_limit(T,df_c)
    #     ca_i_lim=self.cathode_i_limit(T,df_c)        
        
    #     for i in x:
    #         if i < an_i_lim and i< ca_i_lim:
    #             ya.append(self.anodic_over_Kibria_low(T,i)[0])
    #             yca.append(self.cathodic_over_Kibria_low(T,i)[0])
    #             ya_bub.append(self.anodic_over_Bubbles_low(T,i))
    #             yca_bub.append(self.cathodic_over_Bubbles_low(T,i))
                
    #             y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
    #             y_d.append(self.eta_d_Vermeiren(i))
                
    #             tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_low(T,i)[0]+
    #                        self.anodic_over_Bubbles_low(T,i)+
    #                        self.cathodic_over_Kibria_low(T,i)[0]+
    #                        self.cathodic_over_Bubbles_low(T,i)+
    #                        self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i))
    #         elif i < an_i_lim and i> ca_i_lim:
    #             ya.append(self.anodic_over_Kibria_low(T,i)[0])
    #             yca.append(self.cathodic_over_Kibria_high(T,i)[0])
    #             ya_bub.append(self.anodic_over_Bubbles_low(T,i))
    #             yca_bub.append(self.cathodic_over_Bubbles_high(T,i))
    #             y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
    #             y_d.append(self.eta_d_Vermeiren(i))
                
    #             tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_low(T,i)[0]+
    #                        self.anodic_over_Bubbles_low(T,i)+
    #                        self.cathodic_over_Kibria_high(T,i)[0]+
    #                        self.cathodic_over_Bubbles_high(T,i)+
    #                        self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i))
    #         elif i > an_i_lim and i< ca_i_lim:
    #             ya.append(self.anodic_over_Kibria_high(T,i)[0])
    #             yca.append(self.cathodic_over_Kibria_low(T,i)[0])
    #             ya_bub.append(self.anodic_over_Bubbles_high(T,i))
    #             yca_bub.append(self.cathodic_over_Bubbles_low(T,i))
                
    #             y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
    #             y_d.append(self.eta_d_Vermeiren(i))
    #             tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_high(T,i)[0]+
    #                        self.anodic_over_Bubbles_high(T,i)+
    #                        self.cathodic_over_Kibria_low(T,i)[0]+
    #                        self.cathodic_over_Bubbles_low(T,i)+
    #                        self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i))
    #         elif i > an_i_lim and i> ca_i_lim:
    #             ya.append(self.anodic_over_Kibria_high(T,i)[0])
    #             yca.append(self.cathodic_over_Kibria_high(T,i)[0])
    #             ya_bub.append(self.anodic_over_Bubbles_high(T,i))
    #             yca_bub.append(self.cathodic_over_Bubbles_high(T,i))
                
    #             y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
    #             y_d.append(self.eta_d_Vermeiren(i))
    #             tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_high(T,i)[0]+
    #                        self.anodic_over_Bubbles_high(T,i)+
    #                        self.cathodic_over_Kibria_high(T,i)[0]+
    #                        self.cathodic_over_Bubbles_high(T,i)+
    #                        self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i))
    #         else:
    #             pass
                
                
                
                
          
        
    #     df=pd.DataFrame(index=x)
        
       
        
    #     df["Tafel"]=tot
        
       
        
        
    #     return df["Tafel"]
    
    
    
    
    #%%
    def kibria_mix(self,T,COH_conc,delta_l,w,P_s,P):
        x=np.arange(10,2000,10) #mA/cm^2 to be converted in A/cm^2
        ya=[]
        yca=[]
        y_l=[]
        y_d=[]
        tot=[]
        an_i_lim=self.anode_i_limit(T,df_c)
        ca_i_lim=self.cathode_i_limit(T,df_c)        
        
        for i in x:
            if i < an_i_lim and i< ca_i_lim:
                ya.append(self.anodic_over_Kibria_low(T,i)[0])
                yca.append(self.cathodic_over_Kibria_low(T,i)[0])
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_low(T,i)[0]+
                           self.cathodic_over_Kibria_low(T,i)[0]+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i))
            elif i < an_i_lim and i> ca_i_lim:
                ya.append(self.anodic_over_Kibria_low(T,i)[0])
                yca.append(self.cathodic_over_Kibria_high(T,i)[0])
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_low(T,i)[0]+
                           self.cathodic_over_Kibria_high(T,i)[0]+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i))
            elif i > an_i_lim and i< ca_i_lim:
                ya.append(self.anodic_over_Kibria_high(T,i)[0])
                yca.append(self.cathodic_over_Kibria_low(T,i)[0])
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_high(T,i)[0]+
                           self.cathodic_over_Kibria_low(T,i)[0]+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i))
            elif i > an_i_lim and i> ca_i_lim:
                ya.append(self.anodic_over_Kibria_high(T,i)[0])
                yca.append(self.cathodic_over_Kibria_high(T,i)[0])
                y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
                y_d.append(self.eta_d_Vermeiren(i))
                tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_high(T,i)[0]+
                           self.cathodic_over_Kibria_high(T,i)[0]+
                           self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i))
            else:
                pass
                
                
                
                
          
        
        df=pd.DataFrame(index=x)
        
       
        
        df["Tafel"]=tot
        
       
        
        
        return df["Tafel"]
        
     
    
    
    
    #%% Hydroxide concentration evaluation OH abd also o2
    
    def A(self,T_celsius): # temp in Celsius
        A=-0.0032*T_celsius**2-0.1108*T_celsius+1001.7
        return A
    
    
    
    def density(self,A,w):
        rho=A*np.exp(0.86*w)
        return rho
        
    def conc(self,rho,w): #OH concentration 
        concentration=w*rho/56.1 
        return concentration
    
    
    
    #%%
    # O2 concentration
    def f(self,T_kelvin,w,P): #w=0-1 P=[atm] O2 pressure
        M_conc=self.conc(self.density(self.A(T_kelvin-273),w),w)*1000/self.density(self.A(T_kelvin-273),w)
        #M_conc=6.2388
    
        phi=(1/(1+0.102078*M_conc**(1.00044)))**(4.308933)
        f_O2=np.exp((0.046*T_kelvin**2+203.35*T_kelvin*np.log(T_kelvin/298)-(299.378+0.092*T_kelvin)*(T_kelvin-298)-20591)/(8.3144*T_kelvin))
        
        #C_O2=f_O2*phi*P/1000 #mol/dm^3 deviation
    
        C_O2=P*f_O2*phi*self.density(self.A(T_kelvin-273),w)/1000 #mol/dm^3 deviation
    
        return C_O2

#%% reference conditions (MIles 80°C and 50 wt%)
    def ref_conditions(self,T_ref,w_ref,P_ref):

                        
        self.COH_ref=self.conc(self.density(self.A(T_ref-273),w_ref),w_ref)
        self.CO2_ref=self.f(T_ref,w_ref,P_ref)

    
#%% multiple temperature plot
    def multi_T(self,w,COH_conc,P,P_s):
        # plt.figure(5656,figsize=(15,7))
        # ax=plt.gca()
       
       
        kib_mix=pd.DataFrame(index=np.arange(10,2000,10))
        miles_mix=pd.DataFrame(index=np.arange(10,2000,10))
        for i in np.arange(303,373+5,5):
            
            
            kib_mix["temp "+str(i)]=self.kibria_mix(i,COH_conc,self.delta_l,w,P_s,P)
            miles_mix["temp "+str(i)]=self.Miles_mix_bub_act(i,COH_conc,self.delta_l,w,P_s,P)
            # ax.plot(kib_mix.index,kib_mix,label=f"T={i} [K]")

        # ax.legend()
        # ax.set_xlabel("J $[mA/cm^2]$")
        # ax.set_ylabel("$V_{cell} [V]$")
        # ax.grid()
        # ax.set_ylim(1.0, 2.75)
        # ax.set_title("Kibria high Butler-Volmer Temp variation")
   
        return kib_mix,miles_mix
            
    
    
    #%% OH ions mass transport in separator
    
    def eta_sep(self,COH_conc,T,w,P_s,delta_l):
        x=np.linspace(0.5,2000,100) #mA/cm^2 to be converted in A/cm^2
        ya=[]
        yca=[]
        y_l=[]
        y_d=[]
        tot=[]
        T=353 #K
        for i in x:
            ya.append(self.anodic_over_Kibria_high(T,i)[0])
            yca.append(self.cathodic_over_Kibria_high(T,i)[0])
            y_l.append(self.eta_liquid(COH_conc,T,i,delta_l))
            y_d.append(self.eta_d_Vermeiren(i))
            tot.append(self.v_rev(T,w,P_s)+self.anodic_over_Kibria_high(T,i)[0]+self.cathodic_over_Kibria_high(T,i)[0]+self.eta_liquid(COH_conc,T,i,delta_l)+self.eta_d_Vermeiren(i))
        
        v_re=[self.v_rev(T,w,P_s)]*len(x)
        
        df_e=pd.read_csv("Separator voltage experiment.csv",index_col=0)
        
        
        df=pd.DataFrame(index=x)
        #df.index.name("J $[mA/cm^2]$")
        df["$V_{rev}$"]=v_re
        df["$\eta_{an}$"]=[ya[i]+v_re[i] for i in range(len(x))]
        df["$\eta_{ca}$"]=yca
        df["$\eta_l$"]=y_l
        df["$\eta_d$"]=y_d
        # df["$\eta_{an}$"]=[ya[i]+v_re[i] for i in range(len(x))]
        # df["$\eta_{ca}$"]=[yca[i]+ya[i]+v_re[i] for i in range(len(x))]
        # df["$\eta_l$"]=[y_l[i]+yca[i]+ya[i]+v_re[i] for i in range(len(x))]
        # df["$\eta_d$"]=[y_d[i]+y_l[i]+yca[i]+ya[i]+v_re[i] for i in range(len(x))]
        df["$V_{cell}$ or $\eta_d$"]=tot    
        
        
        plt.figure(figsize=(15,7))
        
        plt.plot([i/10 for i in df_e.index],df_e.iloc[:,0].tolist(),"o",c='k',markeredgewidth=3, markersize=12,fillstyle="none", label="Experimental data",linewidth=3)
        plt.plot(df.index,df["$\eta_d$"],marker=".",label="Modelled")
        plt.grid()
        plt.ylim(0,0.34)
        plt.xlabel("J $[mA/cm^2]$")
        plt.xticks(np.arange(0,2000+200,200))
        plt.yticks(np.arange(0,0.34+0.02,0.02))
        plt.legend()
        plt.ylabel("$V_{cell} [V]$")
        plt.tight_layout() 
        plt.title(f"$\eta_d$ with R={self.R_eta_d} $[\Omega \cdot cm^2]$")      

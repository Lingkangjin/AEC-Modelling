# -*- coding: utf-8 -*-
"""
Created on Tue May 24 15:39:53 2022

@author: utente
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pickle



#%%
#solving io with arrenhius plot for tafel





#%% loading parameters

param="Fitting Parameters" #Parameters fodler name

df_c=pd.read_pickle(param+"\df_c.pkl")
DF_ARR=pd.read_pickle(param+"\DF_ARR.pkl")



with open(param+"\\tafel.txt", "rb") as myFile:
    tafel_par = pickle.load(myFile)
    
    

with open(param+"\\alpha.txt", "rb") as myFile:
    alpha_par = pickle.load(myFile)    



with open(param+"\\Vermeiren.txt", "rb") as myFile:
    VE_par = pickle.load(myFile)    

#%% Define the class


class AEC_model(object):
    """
    Alkaline thermo-electrical model, and as main output as the polarisation curve, which different overvoltages are considered:
        -Anodic and Cathodic overvoltage using Tafel semplifications (eta_a and eta_ca) 
        
        -Electrolyte overpotential (eta_l)
        
        -Ion transportation overpotential (eta_d)
        
        -Bubble formation with ohmic behaviour using bubble coverage and Bruggman relationship (eta_bub)
    
    Different sets of conditions of exchange current density (i0) and transfer coefficients (alpha) are provided, based on the literature,
    namely Kibria and Miles references; while for i0 is used Arrenhius fit, transfer coefficients is fitted using a polynomial function,
    in order to get temperature dependency.

    Parameters
    ----------
    T : Float type
        Temperature [K]
        
    COH_conc: float type
        COH concentration [mol/dm^3]
        
    P: float type
        Oxygen pressure [atm]
        
    P_S: float
        System pressure [bar]
        
    delta_l: float  
            reaction distance for ohmic part [cm]

    Returns
    -------
    df: dataframe type
        a dataframe that contains polarisation curve with current density in [mA/cm^2] and V in [V]

    """
    

    def __init__(self, T,P_s,w,delta_l): # initial conditions
        self.T=T       
        self.P_s=P_s
        self.w=w
        self.delta_l=delta_l
    #%%    
    def wt_to_cOH(self,*args):
        
        self.COH_conc=self.conc()
     
    def A(self): 
        """
        function to find the parameter of A, for the empirical expression to calculate concentration of OH        
        
        Returns
        -------
        A : float
            Coefficient which depends on temperature, has the same unit of the denisty [kg/m^3]

        """
        
        
        
        
        A=-0.0032*(self.T-273)**2-0.1108*(self.T-273)+1001.7
        return A
    
    
    
    def density(self):
        """
        Function to calculate the density of the electrolyte based on Gilliam et al. 2007 empricial equation

        Returns
        -------
        rho : float
            Density of the electrolyte which is composed of water and KOH [kg/m^3]

        """
        rho=self.A()*np.exp(0.86*self.w)
        return rho
        
    def conc(self): 
        """
        Calculates the concentration of OH, based on the wt fraction (%) and temperature

        Returns
        -------
        concentration : float
            concentration of the OH in electrolyte [mol/dm^3]

        """
        concentration=self.w*self.density()/56.1  #56.1 is the molecular weight of the KOH [g/mol]
        return concentration    
       
    #%% bubble surface coverage 

    def theta_eps(self,J):
        """
        Function to calculates the bubble effects, 
        in particular the bubble surface coverage (theta) and the diffussivity change (epssilon), with Bruggman relationship

        Parameters
        ----------
        J : float
            Current Density [mA/cm^2]

        Returns
        -------
        theta : float
            bubble surface are coverage
        eps : float
            diffussivity change

        """
        # J_lim=30000   # Piontelli results [mA/cm^2]
        # T_0=273
        #theta=(-97.25+182*(T/T_0)-84*(T/T_0)**2)*(J/J_lim)**(0.3) #To be revised temperature dependence T_0
        theta=0.045*(J)**(0.3) #
        #theta=0.023*(J/J_lim)**(0.3)
        eps=theta*(2/3)
        return theta,eps
    #%% current density limit (changing from low overvoltage to high)
    
    def anode_i_limit(self,df_c):
        """
        Determine the Anode current density limit that the Tafel slope changes, i.e. also io has to be changed

        Parameters
        ----------
        T : float
            Temperature [K]           
        df_c : dataframe
            a dataframe that contains all fitting parameters (2nd order polynomial)

        Returns
        -------
        float
            returns a number which is the current density threshold in [mA/cm^2]

        """
        key="OER"        
        return np.polyval(df_c[key],self.T)
    
    
    def cathode_i_limit(self,df_c):
        """
        Determine the Cathode current density limit that the Tafel slope changes, i.e. also io has to be changed

        Parameters
        ----------
        T : float
            Temperature [K]           
        df_c : dataframe
            a dataframe that contains all fitting parameters (2nd order polynomial)

        Returns
        -------
        float
            returns a number which is the current density threshold in [mA/cm^2]

        """
        key="HER"
        df_c=pd.read_pickle(param+"\\df_c.pkl")    
        return np.polyval(df_c[key],self.T)
    #%% arrenhius exchange current density
    def io_func(self,key):
        """
        Determine the exchangde current density (i0) expressed in [mA/cm^2] using the Arrenhius fit parameters based on the literature experiments

        Parameters
        ----------
        key : String 
            Describes which set of data (Miles or Kibria), which electrode (anode or cathode) and which region (low or high overvoltage region)
        T : float
            Temperature in [K]

        Returns
        -------
        float
            io: exchange current density in [mA/cm^2]

        """
        T_ref=DF_ARR.loc[key][0]
        io_ref=DF_ARR.loc[key][1]
        R=0.00831 #[kJ/mol K]
        E=DF_ARR.loc[key][2]
        
        return io_ref*np.exp((E/R)*((1/T_ref)-(1/self.T)))
   #%% reversible voltage calculation     
    def reve_V(self):
        """
        Empirical equation to determine the reversible part of the voltage

        Parameters
        ----------
        T_in : float
            Temperature [K]

        Returns
        -------
        rev : float
            reversible voltage

        """
        rev=1.5184-1.5421*10**(-3)*self.T+9.523*10**(-5)*self.T*np.log(self.T)+9.84*10**(-8)*self.T**2
        return rev
        
    def v_rev(self):
        """
        Complete reversible overpotential, including the pressure dependence part

        Parameters
        ----------
        T_in : float
            Temperature [K]
        w : float [0-1]
            Weight fraction of the electrolyte 
        P_s : float
            System pressure [bar]

        Returns
        -------
        V_re : float
            Reversible potential

        """
       
        m=self.conc()
        a=-0.0151*m-1.6788*10**(-3)*m**2+2.2588*10**(-5)*m**3
        b=1-1.2062*10**(-3)*m+5.6024*10**(-4)*m**2-7.8228*10**(-6)*m**(3)
        Pv_h20=np.exp(81.6179-7699.68/self.T-10.9*np.log(self.T)+9.5891*10**(-3)*self.T)
        Pv_KOH=np.exp(2.302*a+b*np.log(Pv_h20))
        alfa_H20=np.exp(-0.05192*m+0.003302*m**2+(3.177*m-2.131*m**2)/self.T)
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]
        V_re=self.reve_V()+(R*self.T)/(2*F)*np.log((self.P_s-Pv_KOH)**(1.5)/alfa_H20)
        return V_re  
        
    
        
    #%% anodic overpotential
    
    
    
    def anodic_over_Kibria_low(self,J): #J expressed in [mA/cm^2] 
        """
        Anodic overpotential with Kibria coefficients in low overpotential region

        Parameters
        ----------
        T_in : float
            Temperature [K]-deprecated
        J : float
            Current Density [mA/cm^2]

        Returns
        -------
        eta_an : float
            Anodic overpotential [V]
        i_o : float
            Exchange current density [mA/cm^2]

        """
       
        key= 'Kibria anode-low'
        i_o=self.io_func(key)
        A=np.polyval(tafel_par[key],self.T) 
        eta_an=A*np.log10(J/i_o)
        return eta_an,i_o
    
    def anodic_over_Kibria_high(self,J):#J expressed in [mA/cm^2] 
        """
        Anodic overpotential with Kibria coefficients in high overpotential region

        Parameters
        ----------
        T_in : float
            Temperature [K]-deprecated
        J : float
            Current Density [mA/cm^2]

        Returns
        -------
        eta_an : float
            Anodic overpotential [V]
        i_o : float
            Exchange current density [mA/cm^2]

        """    
    
    
        key= 'Kibria anode-high'
        i_o=self.io_func(key)
        A=np.polyval(tafel_par[key],self.T)         
        eta_an=A*np.log10(J/i_o)
        return eta_an,i_o
    
    def anodic_over_miles_low(self,J):#J expressed in [mA/cm^2]
        """
        Anodic overpotential with Miles coefficients in low overpotential region

        Parameters
        ----------
        T_in : float
            Temperature [K]-deprecated
        J : float
            Current Density [mA/cm^2]

        Returns
        -------
        eta_an : float
            Anodic overpotential [V]
        i_o : float
            Exchange current density [mA/cm^2]

        """
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]
        key='Miles anode-low'
        i_o=self.io_func(key) #Arrenhius
        # i_o=np.polyval(io_par[key],self.T) #polynomial
        alpha=np.polyval(alpha_par[key],self.T)    
        
        eta_an=2.3*R*self.T/(alpha*F)*np.log10(J/i_o)
        return eta_an,i_o    
    
    
    
    
    def anodic_over_miles_high(self,J):#J expressed in [mA/cm^2]
        """
        Anodic overpotential with Miles coefficients in high overpotential region

        Parameters
        ----------
        T_in : float
            Temperature [K]-deprecated
        J : float
            Current Density [mA/cm^2]

        Returns
        -------
        eta_an : float
            Anodic overpotential [V]
        i_o : float
            Exchange current density [mA/cm^2]

        """
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]
        key='Miles anode-high'
        i_o=self.io_func(key)
        
        #i_o=np.polyval(io_par[key],self.T) 
        alpha=np.polyval(alpha_par[key],self.T)
        eta_an=2.3*R*self.T/(alpha*F)*np.log10(J/i_o)
        return eta_an,i_o    
    #%% anodic bubbles overpotential for the activation phase
    
    
    
    def anodic_over_Bubbles_low(self,J):
        """
        Calculate the bubble effect on the activation area on the anode, with Kibria conditions and low overpotential region;
        where since the effective area is less, i.e. the effective current denisty is higher.
        based on the theta, which is the area coverage parameter, assessed with the other function of this class

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an_bu : float
            Anodic overpotential due to bubble formation

        """
        #J expressed in [mA/cm^2]   
    
        key= 'Kibria anode-low'
        i_o=self.io_func(key)
        A=np.polyval(tafel_par[key],self.T) 
        Jeff=J/(1-self.theta_eps(J)[0])
        eta_an_bu=A*np.log10(Jeff/i_o)
        return eta_an_bu
    
    def anodic_over_Bubbles_high(self,J):#J expressed in [mA/cm^2] 
        """
        Calculate the bubble effect on the activation area on the anode, with Kibria conditions and high overvoltage region;
        where since the effective area is less, i.e. the effective current denisty is higher.
        based on the theta, which is the area coverage parameter, assessed with the other function of this class

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an_bu : float
            Anodic overpotential due to bubble formation

        """
        key= 'Kibria anode-high'
        A=np.polyval(tafel_par[key],self.T)     
        i_o=self.io_func(key)
        Jeff=J/(1-self.theta_eps(J)[0])
        eta_an_bu=A*np.log10(Jeff/i_o)
        return eta_an_bu
    
    
    def anodic_over_Bubbles_low_miles(self,J):#J expressed in [mA/cm^2]   
        """
        Calculate the bubble effect on the activation area on the anode, with Miles conditions and low overpotential region;
        where since the effective area is less, i.e. the effective current denisty is higher.
        based on the theta, which is the area coverage parameter, assessed with the other function of this class

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an_bu : float
            Anodic overpotential due to bubble formation

        """
    
        key= 'Miles anode-low'
        i_o=self.io_func(key)
        A=np.polyval(tafel_par[key],self.T) 
        Jeff=J/(1-self.theta_eps(J)[0])
        eta_an_bu=A*np.log10(Jeff/i_o)
        return eta_an_bu
    
    def anodic_over_Bubbles_high_miles(self,J):#J expressed in [mA/cm^2] 
        """
        Calculate the bubble effect on the activation area on the anode, with Miles conditions and high overpotential region;
        where since the effective area is less, i.e. the effective current denisty is higher.
        based on the theta, which is the area coverage parameter, assessed with the other function of this class

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an_bu : float
            Anodic overpotential due to bubble formation

        """
        key= 'Miles anode-high'
        A=np.polyval(tafel_par[key],self.T)     
        i_o=self.io_func(key)
        Jeff=J/(1-self.theta_eps(J)[0])
        eta_an_bu=A*np.log10(Jeff/i_o)
        return eta_an_bu
   
   
    #%% cathodic over potential
    
    def cathodic_over_Kibria_low(self,J):#J expressed in [mA/cm^2]        
        """
        Cathodic overpotential with Kibria coefficients in low overpotential region

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an : float
            Cathodic overpotential due to bubble formation
        i_o : float
            Exchange current density [mA/cm^2]    

        """    
    
        key= 'Kibria cathode-low'
        i_o=self.io_func(key)
        A=np.polyval(tafel_par[key],self.T)    
        eta_an=A*np.log10(J/i_o)
        return eta_an,i_o
    
    def cathodic_over_Kibria_high(self,J):#J expressed in [mA/cm^2]   
        """
        Cathodic overpotential with Kibria coefficients in high overpotential region

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an : float
            Cathodic overpotential due to bubble formation
        i_o : float
            Exchange current density [mA/cm^2]

        """    
        key= 'Kibria cathode-high'
        i_o=self.io_func(key)
        A=np.polyval(tafel_par[key],self.T)    
        eta_an=A*np.log10(J/i_o)
        return eta_an,i_o
    
    def cathodic_over_miles_low(self,J):#J expressed in [mA/cm^2]
        """
        Cathodic overpotential with Miles coefficients in low overpotential region

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an : float
            Cathodic overpotential due to bubble formation
        i_o : float
            Exchange current density [mA/cm^2]

        """    
    
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]
        
        key='Miles cathode-low'
        i_o=self.io_func(key)
        alpha=np.polyval(alpha_par[key],self.T)    
        
        eta_ca=2.3*R*self.T/(alpha*F)*np.log10(J/i_o)
        return eta_ca,i_o
    
    
    def cathodic_over_miles_high(self,J):#J expressed in [mA/cm^2]
        """
        Cathodic overpotential with Miles coefficients in high overpotential region

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an : float
            Cathodic overpotential due to bubble formation
        i_o : float
            Exchange current density [mA/cm^2]

        """    
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]        
        key='Miles cathode-high'
        i_o=self.io_func(key)
        alpha=np.polyval(alpha_par[key],self.T)        
        eta_ca=2.3*R*self.T/(alpha*F)*np.log10(J/i_o)
        return eta_ca,i_o
    #%% cathodic overpotential due to bubbles
    
    def cathodic_over_Bubbles_low(self,J):#J expressed in [mA/cm^2] 
        """
        Calculate the bubble effect on the activation area on the anode, with Kibria conditions and low overpotential region;
        where since the effective area is less, i.e. the effective current denisty is higher.
        based on the theta, which is the area coverage parameter, assessed with the other function of this class

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an_bub : float
            Cathodic overpotential due to bubble formation

        """
        key= 'Kibria cathode-low'
        A=np.polyval(tafel_par[key],self.T)
        i_o=self.io_func(key)
        Jeff=J/(1-self.theta_eps(J)[0])
        eta_an_bub=A*np.log10(Jeff/i_o)
        return eta_an_bub
    
    def cathodic_over_Bubbles_high(self,J):#J expressed in [mA/cm^2]
        """
        Calculate the bubble effect on the activation area on the anode, with Kibria conditions and high overpotential region;
        where since the effective area is less, i.e. the effective current denisty is higher.
        based on the theta, which is the area coverage parameter, assessed with the other function of this class

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an_bub : float
            Anodic overpotential due to bubble formation

        """
        key= 'Kibria cathode-high'
        i_o=self.io_func(key)
        A=np.polyval(tafel_par[key],self.T)    
        Jeff=J/(1-self.theta_eps(J)[0])
        eta_an_bub=A*np.log10(Jeff/i_o)
        return eta_an_bub
    
    def cathodic_over_Bubbles_low_miles(self,J):#J expressed in [mA/cm^2] 
        """
        Calculate the bubble effect on the activation area on the anode, with Miles conditions and low overpotential region;
        where since the effective area is less, i.e. the effective current denisty is higher.
        based on the theta, which is the area coverage parameter, assessed with the other function of this class

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an_bub : float
            Anodic overpotential due to bubble formation

        """
        key= 'Miles cathode-low'
        A=np.polyval(tafel_par[key],self.T)
        i_o=self.io_func(key)
        Jeff=J/(1-self.theta_eps(J)[0])
        eta_an_bub=A*np.log10(Jeff/i_o)
        return eta_an_bub
    
    def cathodic_over_Bubbles_high_miles(self,J):#J expressed in [mA/cm^2]   
        """
        Calculate the bubble effect on the activation area on the anode, with Miles conditions and high overpotential region;
        where since the effective area is less, i.e. the effective current denisty is higher.
        based on the theta, which is the area coverage parameter, assessed with the other function of this class

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_an_bub : float
            Anodic overpotential due to bubble formation

        """
        key= 'Miles cathode-high'
        i_o=self.io_func(key)
        A=np.polyval(tafel_par[key],self.T)    
        Jeff=J/(1-self.theta_eps(J)[0])
        eta_an_bub=A*np.log10(Jeff/i_o)
        return eta_an_bub
 
    #%% Nickel electrode overpotential reference from Master thesis
    
    
    def eta_ni(self,J): # delta_l reaction distance [cm]
        """
        Calculate the Ohmic overpotential due to the ionic conductivity of the electrode, since Nickel is highly conductive, it is neglectable

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta_ni : float
            Electrodes overpotential
        """  
        
        return J*self.delta_l/(1.477*(self.T-273)**2+791.11*(self.T-273)+160285)
    #%% Bubbles ohmic resistance
    
    def eta_ohm_bub(self,J): 
        """
        Evaluation of the ohmic part of the bubble formation, which correlates the surface area coverage (theta) with epsilon, through Bruggman relationship

        Parameters
        ----------
       
        J : float
            Current denisty [mA/cm^2]

        Returns
        -------
        float
            Bubble ohmic overpotential for a given current denisty
       """
        
        
        #con in [mol/dm^3] 
        sigma_l_f=self.sigma_liquid()
        sigma_b_f=sigma_l_f/((1/(1-self.theta_eps(J)[1])**(1.5))-1) #Bruggman relationship [S/cm]
        return J*self.delta_l/(sigma_b_f*1000)
    
    
    
    
    #%% electrolyte ion mass transport overpotential
    
    
    def sigma_liquid(self): 
        """
        Calculates the ionic conductivity of the electrolyte, based on the empirical equation, 
        which depends on the concentration CKOH[mol/dm^3] and the temperature 

        Parameters
        ----------
       

        Returns
        -------
        float
            ionic conductivity [S/cm] or [A/V cm]

        """
        #con=concentration of OH in [mol/dm^3]
        sig_l=-2.041*self.COH_conc-0.0028*self.COH_conc**2+0.005332 *self.COH_conc*self.T+207.2*self.COH_conc/self.T+0.001043*self.COH_conc**3-0.0000003*self.COH_conc**2*self.T**2
        return sig_l #[S/cm] or [A/V*cm]

        
    def eta_liquid(self,J): 
        """
        Electrolyte over voltage calculation for a given current density [mA/cm^2]

        Parameters
        ----------
        J : float 
            Current density [mA/cm^2]

        Returns
        -------
        Electrolyte overvoltage [V]

        """
        
        #J in mA/cm^2
        return(self.delta_l*J/(self.sigma_liquid()*1000)) #1000 is for the conversion mA into A
    
    
   
    #%% ION mass tranport through separator  over voltage eta d

    def eta_d_Vermeiren(self,J): 
        """
        ionic tranport through the diaphram overvoltage;
        with polynomial fit based on Vermeiren experiments, based on 30%wt KOH

        Parameters
        ----------
        J : float
            Current density [mA/cm^2]

        Returns
        -------
        eta : float
           Diaphram overpotential [V]

        """
        
        #VERmeiren model
        R=np.polyval(VE_par["poly"],self.T)[0]  #[ohm * cm^2]
        eta=R*J/1000 #mA/cm^2 into A/cm^2
        return eta
    
    #%% miles high
    
    def Miles(self,x):
        """
        Polarisation curve with Miles condition at high overpotential region.
        Electrode overpotential included
        
        Parameters 
        are all included as initial conditions and arer
        ----------
        T : float 
            Temperature [K]
        COH_conc : float
            OH concentration in electroyte [mol/dm^3], calculated from 
        delta_l : float
            Reaction distance, assumed as sum of the electrode and diaphrame thickness [cm]
        w : float
            Weight fraction of the KOH [0-1]
        P_s : float
            System pressure [bar]
        x : list or array
            List of current denisties [mA/cm^2]

        Returns
        df: dataframe
            dataframe that contains the polarisation curve, where the index are current denisties [mA/cm^2] and column is the tafel caluclated voltage [V]

        """
  
        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        tot=[]
        
        y_ni=[]
       
        
        for i in x:           
            ya.append(self.anodic_over_miles_high(i)[0])
            yca.append(self.cathodic_over_miles_high(i)[0])
            ya_bub.append(self.anodic_over_Bubbles_high_miles(i))
            yca_bub.append(self.cathodic_over_Bubbles_high_miles(i))
            y_ohm_bub.append(self.eta_ohm_bub(i))
            
            y_l.append(self.eta_liquid(i))
            y_d.append(self.eta_d_Vermeiren(i))
            y_ni.append(self.eta_ni(i))
            
            tot.append(self.v_rev()+
                       ya[-1]+
                       yca[-1]+
                       y_l[-1]+
                       y_d[-1]+
                       y_ni[-1])         
    
        df=pd.DataFrame(index=x)
     
        df["Tafel"]=tot
        return df["Tafel"]
   
    
    #%% miles mix
    
    def Miles_bub_act(self,x):
        """
        Polarisation curve with Miles condition at high overpotential region.
        Electrode overpotential included, bubble effects in activation considered
        
        Parameters 
        are all included as initial conditions and arer
        ----------
        T : float 
            Temperature [K]
        COH_conc : float
            OH concentration in electroyte [mol/dm^3], calculated from 
        delta_l : float
            Reaction distance, assumed as sum of the electrode and diaphrame thickness [cm]
        w : float
            Weight fraction of the KOH [0-1]
        P_s : float
            System pressure [bar]
        x : list or array
            List of current denisties [mA/cm^2]

        Returns
        df: dataframe
            dataframe that contains the polarisation curve, where the index are current denisties [mA/cm^2] and column is the tafel caluclated voltage [V]

        """
  
        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        tot=[]
        
        y_ni=[]
       
        
        for i in x:           
            ya.append(self.anodic_over_miles_high(i)[0])
            yca.append(self.cathodic_over_miles_high(i)[0])
            ya_bub.append(self.anodic_over_Bubbles_high_miles(i))
            yca_bub.append(self.cathodic_over_Bubbles_high_miles(i))
            y_ohm_bub.append(self.eta_ohm_bub(i))
            
            y_l.append(self.eta_liquid(i))
            y_d.append(self.eta_d_Vermeiren(i))
            y_ni.append(self.eta_ni(i))
            
            tot.append(self.v_rev()+
                       ya_bub[-1]+
                       yca_bub[-1]+
                       y_l[-1]+
                       y_d[-1]+
                       y_ni[-1])         
    
        df=pd.DataFrame(index=x)
     
        df["Tafel"]=tot
        return df["Tafel"]
    
    #%% Miles ohmic
    
    def Miles_bub_act_ohm(self,x):
        """
        Polarisation curve with Miles condition in high overvoltage region;
        adding the ohmic and activation bubble effects;
        also nichel electrodes overpotential is included

        Parameters
        ----------
        x : list or array
            List of current denisties [mA/cm^2]

        Returns
        df: dataframe
            dataframe that contains the polarisation curve, where the index are current denisties [mA/cm^2] and column is the tafel caluclated voltage [V]

        """
        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        tot=[]
        
        y_ni=[]
        
        for i in x:           
            ya.append(self.anodic_over_miles_high(i)[0])
            yca.append(self.cathodic_over_miles_high(i)[0])
            ya_bub.append(self.anodic_over_Bubbles_high_miles(i))
            yca_bub.append(self.cathodic_over_Bubbles_high_miles(i))
            y_ohm_bub.append(self.eta_ohm_bub(i))
            
            y_l.append(self.eta_liquid(i))
            y_d.append(self.eta_d_Vermeiren(i))
            
            y_ni.append(self.eta_ni(i))
            
            tot.append(self.v_rev()+
                       ya_bub[-1]+
                       yca_bub[-1]+
                       y_l[-1]+
                       y_d[-1]+
                       y_ohm_bub[-1]+
                       y_ni[-1])            
                
                
          
        
        df=pd.DataFrame(index=x)
        
       
        
        df["Tafel"]=tot
        return df["Tafel"]
    
    #%% Kibria high
    def kibria_mix_bub_act(self,x):
        """
        Polarison curve with Kibria conditions and high overpontential region;
        Current denisty threshold applied, two different sets of i0 and alpha deployed
        
        Adding the bubble effects in activation phase

        Parameters
        ----------
        x : list or array
            List of current denisties [mA/cm^2]

        Returns
        -------
        df: dataframe
            dataframe that contains the polarisation curve, where the index are current denisties [mA/cm^2] and column is the tafel caluclated voltage [V]


        """
        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        y_ni=[]
        tot=[]
        
        an_i_lim=self.anode_i_limit(df_c)
        ca_i_lim=self.cathode_i_limit(df_c)
        
        for i in x:
           if i < an_i_lim and i< ca_i_lim:
               ya.append(self.anodic_over_Kibria_low(i)[0])
               yca.append(self.cathodic_over_Kibria_low(i)[0])
               ya_bub.append(self.anodic_over_Bubbles_low(i))
               yca_bub.append(self.cathodic_over_Bubbles_low(i))
               
               y_l.append(self.eta_liquid(i))
               y_d.append(self.eta_d_Vermeiren(i))
               y_ohm_bub.append(self.eta_ohm_bub(i))
               
               y_ni.append(self.eta_ni(i))
               
               
               tot.append(self.v_rev()+
                          ya_bub[-1]+
                          yca_bub[-1]+
                          y_l[-1]+
                          y_d[-1]+
                          # y_ohm_bub[-1]+
                          y_ni[-1])
               
           elif i < an_i_lim and i> ca_i_lim:
               ya.append(self.anodic_over_Kibria_low(i)[0])
               yca.append(self.cathodic_over_Kibria_high(i)[0])
               ya_bub.append(self.anodic_over_Bubbles_low(i))
               yca_bub.append(self.cathodic_over_Bubbles_high(i))
               y_l.append(self.eta_liquid(i))
               y_d.append(self.eta_d_Vermeiren(i))
               y_ohm_bub.append(self.eta_ohm_bub(i))
               
               y_ni.append(self.eta_ni(i))
               
               
               tot.append(self.v_rev()+
                          ya_bub[-1]+
                          yca_bub[-1]+
                          y_l[-1]+
                          y_d[-1]+
                          # y_ohm_bub[-1]+
                          y_ni[-1])
               
           elif i > an_i_lim and i< ca_i_lim:
               ya.append(self.anodic_over_Kibria_high(i)[0])
               yca.append(self.cathodic_over_Kibria_low(i)[0])
               ya_bub.append(self.anodic_over_Bubbles_high(i))
               yca_bub.append(self.cathodic_over_Bubbles_low(i))
               y_ohm_bub.append(self.eta_ohm_bub(i))
               
               y_l.append(self.eta_liquid(i))
               y_d.append(self.eta_d_Vermeiren(i))
               y_ni.append(self.eta_ni(i))
               
               
               tot.append(self.v_rev()+
                          ya_bub[-1]+
                          yca_bub[-1]+
                          y_l[-1]+
                          y_d[-1]+
                          # y_ohm_bub[-1]+
                          y_ni[-1])
               
           elif i > an_i_lim and i> ca_i_lim:
               ya.append(self.anodic_over_Kibria_high(i)[0])
               yca.append(self.cathodic_over_Kibria_high(i)[0])
               ya_bub.append(self.anodic_over_Bubbles_high(i))
               yca_bub.append(self.cathodic_over_Bubbles_high(i))
               y_ohm_bub.append(self.eta_ohm_bub(i))
               
               y_l.append(self.eta_liquid(i))
               y_d.append(self.eta_d_Vermeiren(i))
               
               y_ni.append(self.eta_ni(i))
               
               
               tot.append(self.v_rev()+
                          ya_bub[-1]+
                          yca_bub[-1]+
                          y_l[-1]+
                          y_d[-1]+
                          # y_ohm_bub[-1]+
                          y_ni[-1])
           else:
               pass
       
        
        df=pd.DataFrame(index=x)
        
       
        
        df["Tafel"]=tot
        return df["Tafel"]
    
    
    
    #%% Kibria Mix bubbles activation
    def kibria_mix_bub_act_ohm(self,x):
        """
        Polarison curve with Kibria conditions and high overpontential region;
        different regions considered, where under and upper the defined current density threshold the sets of i0 and alpha change;
        Added the bubble effects in both activation phase and ohmic phase, also electrodes overpotential is included


        Parameters
        ----------
        x : list or array
            List of current denisties [mA/cm^2]

        Returns
        -------
        df: dataframe
            dataframe that contains the polarisation curve, where the index are current denisties [mA/cm^2] and column is the tafel caluclated voltage [V]


        """

        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        tot=[]
        an_i_lim=self.anode_i_limit(df_c)
        ca_i_lim=self.cathode_i_limit(df_c)
        
        y_ni=[]
        
        
        for i in x:
            if i < an_i_lim and i< ca_i_lim:
                ya.append(self.anodic_over_Kibria_low(i)[0])
                yca.append(self.cathodic_over_Kibria_low(i)[0])
                ya_bub.append(self.anodic_over_Bubbles_low(i))
                yca_bub.append(self.cathodic_over_Bubbles_low(i))
                
                y_l.append(self.eta_liquid(i))
                y_d.append(self.eta_d_Vermeiren(i))
                y_ohm_bub.append(self.eta_ohm_bub(i))
                
                y_ni.append(self.eta_ni(i))
                
                
                tot.append(self.v_rev()+
                           ya_bub[-1]+
                           yca_bub[-1]+
                           y_l[-1]+
                           y_d[-1]+
                           y_ohm_bub[-1]+
                           y_ni[-1])
                
            elif i < an_i_lim and i> ca_i_lim:
                ya.append(self.anodic_over_Kibria_low(i)[0])
                yca.append(self.cathodic_over_Kibria_high(i)[0])
                ya_bub.append(self.anodic_over_Bubbles_low(i))
                yca_bub.append(self.cathodic_over_Bubbles_high(i))
                y_l.append(self.eta_liquid(i))
                y_d.append(self.eta_d_Vermeiren(i))
                y_ohm_bub.append(self.eta_ohm_bub(i))
                
                y_ni.append(self.eta_ni(i))
                
                
                tot.append(self.v_rev()+
                           ya_bub[-1]+
                           yca_bub[-1]+
                           y_l[-1]+
                           y_d[-1]+
                           y_ohm_bub[-1]+
                           y_ni[-1])
                
            elif i > an_i_lim and i< ca_i_lim:
                ya.append(self.anodic_over_Kibria_high(i)[0])
                yca.append(self.cathodic_over_Kibria_low(i)[0])
                ya_bub.append(self.anodic_over_Bubbles_high(i))
                yca_bub.append(self.cathodic_over_Bubbles_low(i))
                y_ohm_bub.append(self.eta_ohm_bub(i))
                
                y_l.append(self.eta_liquid(i))
                y_d.append(self.eta_d_Vermeiren(i))
                y_ni.append(self.eta_ni(i))
                
                
                tot.append(self.v_rev()+
                           ya_bub[-1]+
                           yca_bub[-1]+
                           y_l[-1]+
                           y_d[-1]+
                           y_ohm_bub[-1]+
                           y_ni[-1])
                
            elif i > an_i_lim and i> ca_i_lim:
                ya.append(self.anodic_over_Kibria_high(i)[0])
                yca.append(self.cathodic_over_Kibria_high(i)[0])
                ya_bub.append(self.anodic_over_Bubbles_high(i))
                yca_bub.append(self.cathodic_over_Bubbles_high(i))
                y_ohm_bub.append(self.eta_ohm_bub(i))
                
                y_l.append(self.eta_liquid(i))
                y_d.append(self.eta_d_Vermeiren(i))
                
                y_ni.append(self.eta_ni(i))
                
                
                tot.append(self.v_rev()+
                           ya_bub[-1]+
                           yca_bub[-1]+
                           y_l[-1]+
                           y_d[-1]+
                           y_ohm_bub[-1]+
                           y_ni[-1])
            else:
                pass
                
                
                
                
          
        
        df=pd.DataFrame(index=x)
        
       
        
        df["Tafel"]=tot
        return df["Tafel"]
    
    #%% bubbles activation+ ohmic
    
    
    def kibria_bub_act_ohm(self,x):
        """
        Polarisation curve with Kibria conditions and high overvoltage working region; 1 sets of i0 and alphas
        Both activation and ohmic bubble effects are added 
        also the elctrode over potential is added

         Parameters
         ----------
         x : list or array
             List of current denisties [mA/cm^2]

         Returns
         -------
         df: dataframe
             dataframe that contains the polarisation curve, where the index are current denisties [mA/cm^2] and column is the tafel caluclated voltage [V]
        """
        ya=[]
        yca=[]
        ya_bub=[]
        yca_bub=[]
        y_ohm_bub=[]
        
        y_l=[]
        y_d=[]
        tot=[]
        
        y_ni=[]
       
        
        for i in x:
            ya.append(self.anodic_over_Kibria_high(i)[0])
            yca.append(self.cathodic_over_Kibria_high(i)[0])
            ya_bub.append(self.anodic_over_Bubbles_high(i))
            yca_bub.append(self.cathodic_over_Bubbles_high(i))
            y_ohm_bub.append(self.eta_ohm_bub(i))
            
            y_l.append(self.eta_liquid(i))
            y_d.append(self.eta_d_Vermeiren(i))
            
            y_ni.append(self.eta_ni(i))
            
            tot.append(self.v_rev()+                      
                       ya_bub[-1]+  
                       yca_bub[-1]+
                       y_ohm_bub[-1]+
                       y_l[-1]+
                       y_d[-1]+
                       y_ni[-1])
                       
      
        df=pd.DataFrame(index=x)
        
        df["Tafel"]=tot
       
        return df["Tafel"]
    
   
    #%%
    
    def kibria(self,x):
        """
        Polarisation curve with Kibria conditions without bubble effects

        x : list or array
            List of current denisties [mA/cm^2]

        Returns
        -------
        df: dataframe
            dataframe that contains the polarisation curve, where the index are current denisties [mA/cm^2] and column is the tafel caluclated voltage [V]
        """
        ya=[]
        yca=[]
        y_l=[]
        y_d=[]
        tot=[]
    
        
        for i in x:
            ya.append(self.anodic_over_Kibria_high(i)[0])
            yca.append(self.cathodic_over_Kibria_high(i)[0])
            y_l.append(self.eta_liquid(i))
            y_d.append(self.eta_d_Vermeiren(i))
            tot.append(self.v_rev()+
                       ya[-1]+
                       yca[-1]+
                       y_l[-1]+
                       y_d[-1])
            
            
        df=pd.DataFrame(index=x)
        
       
        
        df["Tafel"]=tot
        
       
        
        
        return df["Tafel"]        
    
    

    
    
    
    #%% Not needed anymore (Tafel)
    # # O2 concentration
    # def f(self,T_kelvin,w,P): #w=0-1 P=[atm] O2 pressure
    #     M_conc=self.conc(self.density(self.A(T_kelvin-273),w),w)*1000/self.density(self.A(T_kelvin-273),w)
    #     #M_conc=6.2388
    
    #     phi=(1/(1+0.102078*M_conc**(1.00044)))**(4.308933)
    #     f_O2=np.exp((0.046*T_kelvin**2+203.35*T_kelvin*np.log(T_kelvin/298)-(299.378+0.092*T_kelvin)*(T_kelvin-298)-20591)/(8.3144*T_kelvin))
        
    #     #C_O2=f_O2*phi*P/1000 #mol/dm^3 deviation
    
    #     C_O2=P*f_O2*phi*self.density(self.A(T_kelvin-273),w)/1000 #mol/dm^3 deviation
    
    #     return C_O2


    
  
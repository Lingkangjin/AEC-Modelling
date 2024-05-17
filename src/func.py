# -*- coding: utf-8 -*-
"""
Created on Tue May 24 15:39:53 2022

@author: Lingkang Jin
"""

import numpy as np
import pandas as pd



#%% loading parameters

param="Fitting Parameters" #Parameters folder name



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
    

    def __init__(self, T,P_s,w): # initial conditions
        self.T=T       
        self.P_s=P_s
        self.w=w
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
        theta=((1/30000)**(0.3))*(J)**(0.3) #
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
       
        m=self.conc()*1000/self.density() #update 12/10/22 [mol/kg]
        # m=self.conc()
        a=-0.0151*m-1.6788*10**(-3)*m**2+2.2588*10**(-5)*m**3
        b=1-1.2062*10**(-3)*m+5.6024*10**(-4)*m**2-7.8228*10**(-6)*m**(3)
        Pv_h20=np.exp(81.6179-7699.68/self.T-10.9*np.log(self.T)+9.5891*10**(-3)*self.T)
        Pv_KOH=np.exp(2.302*a+b*np.log(Pv_h20))
        alfa_H20=np.exp(-0.05192*m+0.003302*m**2+(3.177*m-2.131*m**2)/self.T)
        R=8.31 # [J/mol*K]
        F=96500 #[C/mol]
        V_re=self.reve_V()+(R*self.T)/(2*F)*np.log((self.P_s-Pv_KOH)**(1.5)/alfa_H20)
        return V_re  
        
    
        
    
 
   
    
    
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
    
    
   
   
    
  
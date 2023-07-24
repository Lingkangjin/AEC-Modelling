# -*- coding: utf-8 -*-
"""
Created on Mon May 30 11:04:59 2022

@author: utente
"""

import numpy as np
#%%

#models valid if temperature less than 1000 K

class Thermal_properties(object): #Shomate
    """Compute thermal properties of substances based on empirical formulations."""


    def __init__(self, T_kel):
        """Initialize the Thermal_properties object.

       Args:
           T_kel (float): Temperature in Kelvin.

       Returns:
           None
       """
        self.T_kel=T_kel
    
   
        
    def water(self,*args):
        """Calculate thermal properties of water.

        Args:
            *args: Additional arguments (not used in the method).

        Returns:
            tuple: Tuple containing the specific heat capacity (cp), enthalpy (h), and entropy (s).
        """
        T=self.T_kel/1000
        A=-203.6
        B=1523.29
        C=-3196.413
        D=2474.455
        E=3.855326
        F=-256.5478
        G=-488.7163
        H=-285.8304
        h=round(A*T+B*T**2/2+C*T**3/3+D*T**4/4-E/T+F,2)*1000 #J/mol #deleted h
        cp=A+B*T+C*T**2+D*T**3+E/(T**2) #J/mol K      
        s=A*np.log(T)+B*T+C*T**2/2+D*T**3/3-E/(2*T**2)+G
        return cp,h,s

    def hydrogen(self,*args): #j/K*mol Shomate equations
        """Calculate thermal properties of hydrogen.

        Args:
            *args: Additional arguments (not used in the method).

        Returns:
            tuple: Tuple containing the specific heat capacity (cp), total enthalpy (h_tot), and entropy (s).
        """
        T=self.T_kel/1000
        A=33.066178
        B=-11.363417
        C=11.432816
        D=-2.772874
        E=-0.158558
        F=-9.980797
        G=172.707974
        H=0
        cp=A+B*T+C*T**2+D*T**3+E/(T**2) #J/mol K
        h=round(A*T+B*T**2/2+C*T**3/3+D*T**4/4-E/T+F-H,2)*1000 #J/mol
        
        h_h2_298=0  #J/mol
        s=A*np.log(T)+B*T+C*T**2/2+D*T**3/3-E/(2*T**2)+G
        
        # h_h2_298=8468  #J/mol
        h_tot=h+h_h2_298
        return cp,h_tot,s

    def oxygen(self,*args):
        """Calculate thermal properties of oxygen.

        Args:
            *args: Additional arguments (not used in the method).

        Returns:
            tuple: Tuple containing the specific heat capacity (cp), total enthalpy (h_tot), and entropy (s).
        """
        T=self.T_kel/1000
        A=31.32234
        B=-20.23531
        C=57.86644
        D=-36.50624
        E=-0.007374
        F=-8.903471
        G=246.7945
        H=0
        cp=A+B*T+C*T**2+D*T**3+E/(T**2) #J/mol K
        h=round(A*T+B*T**2/2+C*T**3/3+D*T**4/4-E/T+F-H,2)*1000 #J/mol
        # h_O2_298=8680 #J/mol
        h_O2_298=0
        h_tot=h+h_O2_298
        s=A*np.log(T)+B*T+C*T**2/2+D*T**3/3-E/(2*T**2)+G
        return cp,h_tot,s    
    
    
    
    def delta_h(self):
        """Calculate the change in enthalpy for the electrolysis of water.
    
          Returns:
              float: Change in enthalpy.
        """
        return(-self.water(self.T_kel)[1]+self.hydrogen(self.T_kel)[1]+0.5*self.oxygen(self.T_kel)[1])
    
    
    def delta_s(self):
        """Calculate the change in entropy for the electrolysis of water.

       Returns:
           float: Change in entropy.
       """
        return(-self.water(self.T_kel)[2]+self.hydrogen(self.T_kel)[2]+0.5*self.oxygen(self.T_kel)[2])
import numpy as np
import CoolProp.CoolProp as CP

def void_fraction(self,P_r,x_r):
    '''
        Zivi model (1963), successful in predicting DP (Bae et al. 1969) 
        and HTC (Bae et al. 1969; Soliman et al. 1968) during condensation
    
    '''
    rho_f=CP.PropsSI('D','P',P_r*1e3,'Q',0,self.fld)        # [kg/m^3]
    rho_g=CP.PropsSI('D','P',P_r*1e3,'Q',1,self.fld)
    S = (rho_f / rho_g)**(1/3)
    alpha = 1 / (1 + ((1 - x_r) / x_r) * (rho_g / rho_f) * S)
    rho_r = rho_g*alpha+rho_f*(1-alpha)

    return  alpha,rho_r
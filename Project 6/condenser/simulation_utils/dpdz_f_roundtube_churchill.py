import numpy as np
import CoolProp.CoolProp as CP

def dpdz_f_roundtube_churchill(Re,Dh,roughness):
    '''
        Calculate friction factor of H2OEG
        Churchill (1977) correlation
         * Input Variables:
            @Re         : Reynolds number
            @Dh         : hydraulic diameter
            @sigma      : viscosity
        
         * Output variable:
            @dpdz_f     : [Pa]
    '''

    A = (-2.456*np.log((7/Re)**0.9+0.27*roughness/Dh))**16; B = (37530/Re)**16
    dpdz_f = 8*((8/Re)**12+1/(A+B)**1.5)**(1/12)

    return dpdz_f
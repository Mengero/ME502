import numpy as np
from scipy import interpolate
import CoolProp.CoolProp as CP

def htc_H2OEG(Re,T_r,P_r,W_port,H_port,Dh_port):
    '''
        !Ciatation: from lecture
         * Input Variables:
            @Re         : Reynolds number
            @T_r        : refrigerant temperature, [C]
            @P_r        : refriegrant pressure, [kPa]
            @W_port     : port width, [m]
            @H_port     : port height, [m]
            @Dh_port    : hydraulic diameter of port, [m]
        
         * Output Variables:
            @htc        : heat transfer coefficient of H2OEG, [W/m^2-K]
    '''
    fld = 'INCOMP::MEG[0.5]'            # change the ratio if needed, currently 50%
    
    if Re<2300:     # laminar for roundtube
        # interpolation data from table
        ratio_ab = np.array([1.0,1.43,2.0,3.0,4.0,8.0])
        Nu_D_q = np.array([3.61,3.73,4.12,4.79,5.33,6.49])
        Nu_D_Ts = np.array([2.98,3.08,3.39,3.96,4.44,5.60])
        Nu_D_ave = (Nu_D_q+Nu_D_Ts)/2
        htc_H2OEG_Laminar=interpolate.UnivariateSpline(ratio_ab,Nu_D_ave,s=0)        # pass all the points forceabley
        
        Nu_D = htc_H2OEG_Laminar(W_port/H_port)
    else:           # turbulent flow
        # Gnielinski correlation
        f = (0.79*np.log(Re)-1.64)**(-2)
        Pr = CP.PropsSI('PRANDTL','T',T_r+273.15,'P',P_r*1e3,fld)
        Nu_D = (f/8)*(Re-1000)*Pr/(1+12.7*(f/8)**(1/2)*(Pr**(2/3)-1))

    k = CP.PropsSI('L','T',T_r+273.15,'P',P_r*1e3,fld)     # [W/m-K]
    htc = Nu_D*k/Dh_port

    return htc
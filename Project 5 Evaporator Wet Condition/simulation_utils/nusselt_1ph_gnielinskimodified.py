import numpy as np
from scipy import interpolate
import CoolProp.CoolProp as CP

def nusselt_1ph_gnielinskimodified(Re_D,Pr,AspRat):
    '''
        !Citation: Gnielinsk, 1976 Correlation for heat transfer inside tubes
        !Citation: Kandlikar et al, 2014, Heat Transfer and Fluid Flow in Minichannels and Microchannels, 2nd Edition, Chapter 3.
        Kandlikar and Steinke 2003 recommended to use fully developed laminar for Re_D <= 1600, Gnielinski for Re_D >= 3000, and interpolation for 1600 < Re_D < 3000

         * Input Variables
            @Re_D                   : Reynolds number based on hydraulic diameter
            @Pr                     : Prandtl number
            @AspRat                 : Cross-section aspect ratio, 0 for circulat, 0~1 for rectangular
         * Output Variables:
            @Nu_D                   : Nusselt number based on hydraulic diameter
    '''
    if Re_D<=2000:          # Laminar fully-developed, constant heat flux
        if AspRat==0:
            Nu_D = 4.364
        else:
            if AspRat > 1:
                AspRat = 1/AspRat
            ratio_ab = np.array([1.0,1.43,2.0,3.0,4.0,8.0])
            Nu_D_q = np.array([3.61,3.73,4.12,4.79,5.33,6.49])
            Nu_D_Ts = np.array([2.98,3.08,3.39,3.96,4.44,5.60])
            Nu_D_ave = (Nu_D_q+Nu_D_Ts)/2
            Nu_D_Laminar=interpolate.UnivariateSpline(ratio_ab,Nu_D_ave,s=0)        # pass all the points forceabley
            Nu_D = Nu_D_Laminar(AspRat)
    elif Re_D>2200 and Re_D<2400:           # Transition from laminar to turbulent
        Nu_D = 0.08436*Re_D + 408302/Re_D - 366.82              # Interpolation for the transition region, from Jun's model
    else:           # Re_D >= 2400, fully turbulent
        f = (0.79*np.log(Re_D)-1.64)**(-2)
        Nu_D = (f/8)*(Re_D-1000)*Pr/(1+12.7*(f/8)**(1/2)*(Pr**(2/3)-1))
        
    return Nu_D
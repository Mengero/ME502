import numpy as np
import CoolProp.CoolProp as CP
from simulation_utils.nusselt_1ph_gnielinskimodified import nusselt_1ph_gnielinskimodified

def htc_2ph_Kim_2013(Ref,q_H,G,T_sat,x,Dh,P_H,P_F,AspRat):
    '''
        !Citation: Kim, Sung-Min, and Issam Mudawar. "Universal approach to predicting saturated flow boiling heat transfer in mini/micro-channels–Part II. Two-phase heat transfer coefficient." International Journal of Heat and Mass Transfer 64 (2013): 1239-1256.
         * Input variables
            @q_H                    : effective heat flux, [W/m^2] (averaged heat flux)
            @G                      : refrigerant mass flux, [kg/s]
            @Dh                     : hydraulic diameter, [m]
            @P_H                    : heated perimeter of channel, [m] (outer perimeter)
            @P_F                    : wetted perimeter of channel, [m] (inner perimeter)

         * Output variables
            @htc_tp                 : [W/m^2-K]
    '''
    
    P = CP.PropsSI('P','T',T_sat+273.15,'Q',x,Ref)/1e3      # [kPa]
    P_crit = CP.PropsSI('PCRIT','T',T_sat+273.15,'Q',x,Ref)/1e3     # [kPa]
    h_f = CP.PropsSI('H','T',T_sat+273.15,'Q',0,Ref)        # [J/kg]
    h_g = CP.PropsSI('H','T',T_sat+273.15,'Q',1,Ref)
    h_fg = h_g - h_f
    mu_f=CP.PropsSI('V','T',T_sat+273.15,'Q',0,Ref)         # [N/m^2-s]
    mu_g=CP.PropsSI('V','T',T_sat+273.15,'Q',1,Ref)
    sigma=CP.PropsSI('I','T',T_sat+273.15,'Q',x,Ref)
    rho_f=CP.PropsSI('D','T',T_sat+273.15,'Q',0,Ref)
    rho_g=CP.PropsSI('D','T',T_sat+273.15,'Q',0,Ref)
    Pr_f=CP.PropsSI('PRANDTL','T',T_sat+273.15,'Q',0,Ref)
    k_f=CP.PropsSI('L','T',T_sat+273.15,'Q',0,Ref)


    Bo = q_H/(G*h_fg)
    P_R = P/P_crit
    Re_f = G*(1-x)*Dh/mu_f
    We_fo = G**2*Dh/(rho_f*sigma)
    X_tt = (mu_f/mu_g)**0.1 * ((1-x)/x)**0.9 * (rho_g/rho_f)**0.5
    Ca = mu_f * G / (rho_f * sigma)

    x_di = 1.4 * We_fo**0.03 * P_R**0.08 - 15.0 * (Bo*P_H/P_F)**0.15 * Ca**0.35 * (rho_g/rho_f)**0.06

    x_tp = min(x,x_di)

    Re_f = G*(1-x_tp)*Dh/mu_f
    X_tt = (mu_f/mu_g)**0.1 * ((1-x_tp)/x_tp)**0.9 * (rho_g/rho_f)**0.5
    htc_cb = (5.2 * (Bo*P_H/P_F)**0.08 * We_fo**(-0.54) + 3.5 * (1/X_tt)**0.94 * (rho_g/rho_f)**0.25) * (0.023*Re_f**0.8 * Pr_f**0.4 * k_f/Dh)
    htc_nb = (2345 * ( Bo*P_H/P_F)**0.7 * P_R**0.38 * (1-x)**(-0.51)) * (0.023 * Re_f**0.8 * Pr_f**0.4 * k_f/Dh)
    htc_tp = np.sqrt(htc_nb**2 + htc_cb**2)

    if x > x_di:
        Re_g = G*Dh/mu_g
        Pr_g = CP.PropsSI('PRANDTL','T',T_sat+273.15,'Q',1,Ref)
        Nu_g = nusselt_1ph_gnielinskimodified(Re_g,Pr_g,AspRat)
        k_g = CP.PropsSI('L','T',T_sat+273.15,'Q',1,Ref)      # [W/m-K]
        htc_g = Nu_g * k_g / Dh
        htc_tp = (htc_tp * (1 - x) + htc_g * (x - x_di)) / (1-x_di)

    return htc_tp

# q_H = 2.4*1e4           # [W/m^2]
# Ref = 'R1234yf'
# Dh=0.002       # [m]
# T_sat = 15      # [C]
# x = 0.31
# P_H = 0.02194070751110265 
# P_F = 0.027541592653589794
# G = 200     # [kg/m^3]
# print(htc_2ph_Kim_2013(q_H,G,T_sat,x,Dh,P_H,P_F))
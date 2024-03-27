import numpy as np
import CoolProp.CoolProp as CP

def htc_2ph_y1234yf_shah(fld,G_ref,Dh_port,P_r,x_r):

    '''
        !Citation: Shah, Mirza M. "Improved correlation for heat transfer during condensation in conventional and mini/micro channels." International Journal of Refrigeration 98 (2019): 222-237.
         * Input Variables:
            @P_r            : refrigerant pressure, [kPa]
            @x_r            : refrigerant vapor quality, [0-1]

         * Output Variables:
            @htc_ref        : refrigerant heat transfer coefficient, [W/m^2-K]
    '''

    mu_l = CP.PropsSI('V','P',P_r*1e3,'Q',0,fld)       # [N/m^2-s]
    mu_g = CP.PropsSI('V','P',P_r*1e3,'Q',1,fld)
    rho_l = CP.PropsSI('D','P',P_r*1e3,'Q',0,fld)      # [kg/m^3]
    rho_g = CP.PropsSI('D','P',P_r*1e3,'Q',1,fld)
    Pr_l = CP.PropsSI('PRANDTL','P',P_r*1e3,'Q',0,fld)
    k_l = CP.PropsSI('L','P',P_r*1e3,'Q',0,fld)        # [W/m-K]
    
    Re_LT = G_ref*Dh_port/mu_l
    h_LT = 0.023*Re_LT**0.8*Pr_l**0.4*k_l/Dh_port
    htc_ref = h_LT*(1 + 1.128* x_r**0.817 * (rho_l/rho_g)**(0.3685) * (mu_l/mu_g)**0.2363 * (1-mu_g/mu_l)**2.144 * Pr_l**(-0.1))
    
    return htc_ref
import numpy as np
import CoolProp.CoolProp as CP
from simulation_utils.nusselt_1ph_gnielinskimodified import nusselt_1ph_gnielinskimodified

def htc_2ph_chen_1966(Ref,T_r,P_r,x_r,G,Dh,R_air,T_ai,A_ref,AspRat):
    '''
        !Citation: S. G. Kandlikar, 1990. A general correlation for saturated two-phase flow boiling heat transfer inside horizontal and vertical tubes. ASME J. of Heat Transfer 112:219-228
        !Citation: S. G. Kandlikar, 1990. Flow boiling maps for water, R-22, and R-134a in the saturated region, Int. Heat Transfer Conference, Jurisalem
        !Citation: S. G. Kandlikar, Mark E. Steinke, Predicting heat transfer during flow boiling in minichannels and microcchannels. ASHRAE Transactions CH-03-13-1, 2003

         * Input Variables
            @Ref            : name of refrigerant
            @T_r            : refrigerant temperature, [C]
            @P_r            : refrigerant pressure, [kPa]
            @x_r            : refrigerant quality, [-]
            @G              : refrigerant mass flux, [kg/m**2-s]
            @Dh             : hydraulic diameter, [m]
            @R_air          : air side thermal resistance, [K/W]
            @T_ai           : air side temperature, [C]
            @A_ref          : refrigerant side heat transfer area, [m**2]
            @AspRat         : aspect ratio, [-]

         * Output Variables
            @htc_tp         : two phase heat tranfer coefficient, [W/m**2-K]
            @T_wall         : tube wall temperature, [C]
    '''

    x_r = min(max(x_r,1e-4),0.9999)         # avoid singularity
    rho_f = CP.PropsSI('D','P',P_r*1e3,'Q',0,Ref)               # [kg/m**3]
    rho_g = CP.PropsSI('D','P',P_r*1e3,'Q',1,Ref)
    mu_f = CP.PropsSI('V','P',P_r*1e3,'Q',0,Ref)                # [N/m**2-s]
    mu_g = CP.PropsSI('V','P',P_r*1e3,'Q',1,Ref)
    k_f = CP.PropsSI('L','P',P_r*1e3,'Q',0,Ref)                 # [W/m-K]
    sigma = CP.PropsSI('I','P',P_r*1e3,'Q',x_r,Ref)             # [N/m]
    Cp_f = CP.PropsSI('C','P',P_r*1e3,'Q',0,Ref)                # [J/kg-K]
    h_f = CP.PropsSI('H','P',P_r*1e3,'Q',0,Ref)                 # [J/kg]
    h_g = CP.PropsSI('H','P',P_r*1e3,'Q',1,Ref)
    h_fg = h_g - h_f
    Pr_f = CP.PropsSI('PRANDTL','P',P_r*1e3,'Q',0,Ref)          # [-]

    x_r_tp = x_r
    x_r = min(x_r,0.7)
    
    Re_f = G * (1 - x_r) * Dh / mu_f
    X_tt = ((1 - x_r) / x_r)**0.9 * (rho_g / rho_f)**0.5 * (mu_f / mu_g)**0.1      # Zhang et al 2004 Eq (7)
    if X_tt >= 10:
        F=1     # Zhang et al 2004 Eq (3)
    else:
        F=2.35 * (1 / X_tt + 0.213)**0.736           # Zhang et al 2004 Eq (2)

    htc_sp = 0.023 * Re_f**0.8 * Pr_f**0.4 * (k_f / Dh)   # Zhang et al 2004 Eq (4), Dittus-Boelter equation
    Re_tp = F**1.25*Re_f
    S = 1 / (1 + 2.53e-6 * Re_tp**1.17)
    T_wall = (T_r + T_ai) / 2
    T_wall_0 = T_wall+1
    while abs(T_wall - T_wall_0) >= 0.01:
        T_wall_0 = T_wall
        DELTAT_sat = T_wall - T_r
        DELTAP_sat = CP.PropsSI('P','T',T_wall+273.15,'Q',0.5,Ref)/1e3-CP.PropsSI('P','T',T_r+273.15,'Q',0.5,Ref)/1e3
        htc_nb = 0.00122 * (k_f**0.79 * Cp_f**0.45 * rho_f**0.49 / (sigma**0.5 * mu_f**0.29 * h_fg**0.24 * rho_g**0.24)) * DELTAT_sat**0.24 * DELTAP_sat**0.75    # Zhang et al 2004 Eq (6), Forster and Zuber equation
        htc_tp = S * htc_nb + F * htc_sp		# Zhang et al 2004 Eq (1)
        T_wall = (1 / R_air * T_ai + htc_tp * A_ref * T_r) / (1 / R_air + htc_tp * A_ref)
    
    if x_r_tp>0.7:           # for x_r > 0.7, assume it is dry out region (high quality region), use linear interpolation between single phase htc and 2 phase htc
        Re_g = G*Dh/mu_g
        Pr_g = CP.PropsSI('PRANDTL','P',P_r*1e3,'Q',1,Ref)
        Nu_g = nusselt_1ph_gnielinskimodified(Re_g,Pr_g,AspRat)
        k_g = CP.PropsSI('L','T',T_r+273.15,'Q',1,Ref)      # [W/m-K]
        htc_g = Nu_g * k_g / Dh
        htc_tp = (htc_tp * (1 - x_r_tp) + htc_g * (x_r_tp - 0.7)) / 0.3
        T_wall = (1 / R_air * T_ai + htc_tp * A_ref * T_r) / (1 / R_air + htc_tp * A_ref)

    return htc_tp, T_wall
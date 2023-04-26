import CoolProp.CoolProp as CP
import numpy as np

'''
    Assumptions:
        1. Uniform air-side surface temperature, T_sa
        2. Refrigerant side surface temperature is uniform, T_sr
        3. Condensation happens when T_sa lower than T_dp
        4. Water film is saturated water at T_sa
        5. Humid air at water film is saturated air, T_sa, and RH=1
        6. Convection heat transfer htc_a, htc_r are constant
'''

def dehumidification(T_ai,T_ri,T_sa,htc_a,A_sa,cp_a,rho_a,Le_a,rho_As,\
                     Vol_dot_a,Vda_i,m_dot_da,omega_Ai,h_fg,m_dot_r,cp_r,R_wall,\
                     htc_r,A_sr,rho_Ai,eta_t):
    '''
     * Input Variables:
        @T_ai               : air inlet temperature, [C]
        @T_ri               : refrigerant inlet temperature, [C]
        @T_sa               : air side surface temperature, [C]
        @htc_a              : air side heat transfer coefficient, [kW/m^2-K]
        @A_sa               : air side surface area, [m^2]
        @cp_a               : air heat capacity, [kJ/kg-K]
        @rho_ai             : air inlet density, [kg/m^3]
        @Le_a               : Lewis number, Le_a = alpha_a/D_AB_a, D_AB_a = 0.26e-4
        @rho_as             : air density at surface temperature, [kg/m^3]
        @Vol_dot_a          : volumetric flow rate of air, [m^3/s]
        @Vda_i              : inlet dry air specific volume, [m^3/kg]
        @m_dot_da           : dry air mass flow rate, [kg/s]
        @omega_ai           : air inlet humidity ratio, [-]
        @h_fg               : latent heat of liquid-vapor phase change of water at T_sa
        @m_dot_r            : refrigerant mass flow rate, [kg/s]
        @cp_r               : specific heat of refrigerant, [kJ/kg-K]
        @R_wall             : tube wall thermal resistance
        @htc_r              : refrigerant heat transfer coefficient, [kW/m^2-K]
        @A_sr               : refrigerant side heat transfer area, [m^2]
        @eta_t              : total surface efficiency, [-]

     * Output Variables:
        @Q_a_sen            : air side sensible heat transfer, [kW]
        @Q_a_lat            : air side latent heat transfer, [kW]
        @Q_r                : refrigerant side heat transfer, [kW]
    
    '''
    # Q_a_lat+Q_a_sen = Q_r

    # ---------------- Q_a_sen Calculation ----------------- #
    T_ao = T_sa + (T_ai-T_sa) * np.exp( (-htc_a*A_sa*eta_t) / (m_dot_da*cp_a))
    Q_a_sen = m_dot_da*cp_a*(T_ai-T_ao)
    
    # ---------------- Q_a_lat Calculation ----------------- #
    n = 1/3         # from heat and mass transfer analogy
    mtc_a = htc_a / (rho_a * cp_a * Le_a**(1-n))            # [m/s]
    rho_Ao = rho_As + (rho_Ai - rho_As) * np.exp(-mtc_a*A_sa*eta_t**0.5/Vol_dot_a)
    omega_Ao = rho_Ao*Vda_i
    m_dot_cond = m_dot_da * (omega_Ai - omega_Ao)
    Q_a_lat = m_dot_cond*h_fg

    # ---------------- Q_r Calculation ----------------- #
    T_ro = T_sa - (T_sa - T_ri) * np.exp(-1 / (m_dot_r * cp_r * (R_wall + 1 / (htc_r * A_sr))))
    Q_r = m_dot_r * cp_r * (T_ro - T_ri)

    return Q_a_sen, Q_a_lat, Q_r, omega_Ao
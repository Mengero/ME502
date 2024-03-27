import numpy as np
import CoolProp.CoolProp as CP
import simulation_utils.DP_ref as DP_ref
from simulation_utils.nusselt_1ph_gnielinskimodified import nusselt_1ph_gnielinskimodified
from dehumidification_model import dehumidification
from simulation_utils.htc_2ph_Kim_2013 import htc_2ph_Kim_2013
from simulation_utils.DP_air import DP_air

def epsilon_NTU(self,x_ro,h_ro,P_ro,m_dot_ref,T_ai,P_ai,RH_ai,m_dot_air,m_dot_da,R_air,htc_air,eta_t,G_ref_elem,N_slab):
    '''
        epsilon-NTU method
        Since this is for cross flow and counter flow evaporator, we need to back iterate the whole process

         * Input Variables
            @x_ro           : refrigerant outlet vapor quality, [-]
            @h_ro           : refrigerant outlet enthalpy, [J/kg]
            @P_ro           : refrigerant outlet pressure, [kPa]
            @m_dot_ref      : refrigerant mass flow rate, [kg/s]
            @T_ai           : air inlet temperature, [C]
            @P_ai           : air inlet pressure, [kPa]
            @RH_ai          : air inlet relative humidity, [-]
            @m_dot_air      : air mass flow rate, [kg/s]
            @m_dot_da       : dry air mass flow rate
            @R_air          : air thermal resistance, [K/kW-m^2]
            @htc_air        : air side heat transfer coefficient, [W/m^2-K]
            @eta_t          : surface total efficiency (energy), [-]
            @G_ref_elem     : refrigerant mass velocity, [kg/m^2-s]
            @htc_correlation    : htc correlation method.
        
         * Output Variables
            @NTU, sigma
            @Q              : heat transfer rate, [W]
            @T_ro           : refrigerant outlet temperature, [C]
            @T_ao           : air outlet temperature, [C]
    '''
    if N_slab%2 == 0:
        z_i = 0; z_o = self.L_elem
    else:
        z_i = self.L_elem; z_o = 0

    # Refrigerant Properties
    k_ref = CP.PropsSI('L','H',h_ro,'P',P_ro*1e3,self.fld)              # [W/m-K]
    T_ro = CP.PropsSI('T','P',P_ro*1e3,'H',h_ro,self.fld)-273.15        # [C]
    cp_ref = CP.PropsSI('C','H',h_ro,'P',P_ro*1e3,self.fld)             # [J/kg-K]
    rho_ref=CP.PropsSI('D','H',h_ro,'P',P_ro*1e3,self.fld)
    C_ref = cp_ref*m_dot_ref            # [J/K]
    Pr = CP.PropsSI('PRANDTL','H',h_ro,'P',P_ro*1e3,self.fld)        # [-]
    mu_ref = CP.PropsSI('V','H',h_ro,'P',P_ro*1e3,self.fld)              # [N/m^2-s]

    # Iteration with dehumidification model
    Q_r = 9999
    Q_a = 0
    Tdp_ai = CP.HAPropsSI('Tdp','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)-273.15
    T_sa = Tdp_ai-1

    # Air side Properties
    cp_air = CP.HAPropsSI('C','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)
    C_air = cp_air*m_dot_da
    cp_ha = CP.HAPropsSI('Cha','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)
    C_ha = cp_ha*m_dot_air
    rho_air = 1/CP.HAPropsSI('Vha','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)          # [kg/m^3]
    k_air = CP.HAPropsSI('K','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)      # [W/m-K]
    alpha_a = k_air / (rho_air*cp_air)               # [m^2/s]
    D_AB_a = 0.26e-4                                # [m^2/s]
    Le_a = alpha_a / D_AB_a
    Vda_i = CP.HAPropsSI('Vda','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)         # [m^3/kg]
    omega_Ai = CP.HAPropsSI('W','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)        # [-]
    rho_Ai = omega_Ai / Vda_i
    
    if (x_ro>1) or (x_ro<0):          # SP
        Re_ref = self.G_ref_elem*self.Dh_port/mu_ref
        Nu_1ph = nusselt_1ph_gnielinskimodified(Re_ref,Pr,self.asp_ratio)
        htc_ref = Nu_1ph*k_ref/self.Dh_port            # [W/m^2-K]
        C_min = min(C_ref,C_ha)
        C_max = max(C_ref,C_ha)
        C_ratio = C_min/C_max
        while abs(Q_r-Q_a)>1e-3:
            h_f = CP.PropsSI('H','T',T_sa+273.15,'Q',0,'H2O')      # [J/kg]
            h_g = CP.PropsSI('H','T',T_sa+273.15,'Q',1,'H2O')
            h_fg = h_g - h_f
            if T_sa < Tdp_ai:               # apply 0.1 is to avoid singularity, or else the water will be evaporated to air flow which makes no sense
                omega_s = CP.HAPropsSI('W','T',T_sa+273.15,'P',P_ai*1e3,'RH',1)             # [-]
                Vda_s = CP.HAPropsSI('Vda','T',T_sa+273.15,'P',P_ai*1e3,'RH',1)             # [m^3/kg]
                rho_As = omega_s / Vda_s            # [kg/m^3], kg of H2O/m^3
                [Q_a_sen, Q_a_lat, Q_r, omega_Ao, T_ao] = dehumidification(T_ai,T_ro,x_ro,T_sa,htc_air,self.A_air_elem,cp_air,rho_air,Le_a,rho_As,\
                                                                     self.vol_dot_air_elem,Vda_i,m_dot_da,omega_Ai,h_fg,self.m_dot_ref_elem,cp_ref,self.R_wall,\
                                                                     htc_ref,self.A_ref_elem,rho_Ai,eta_t)
                Q_a = Q_a_sen+Q_a_lat
                Q = Q_a
            else:
                omega_Ao = omega_Ai
                UA = 1 / (R_air + 1/(htc_ref*self.A_ref_elem))
                NTU = UA/C_min
                sigma = 1-np.exp((1/C_ratio)*NTU**0.22*(np.exp(-C_ratio*NTU**0.78)-1))

                T_ri = (C_ref*T_ro-sigma*C_min*T_ai)/(C_ref - sigma*C_min)
                Q = C_ref*(T_ro-T_ri)
                Q_a = Q; Q_r = Q
                Q_a_sen = Q_a; Q_a_lat = 0
                T_ao = T_ai - Q/C_ha
            
            T_sa += (Q_a-Q_r)/10
        h_ri = h_ro - Q/m_dot_ref
        Vel_i=self.G_ref_elem/rho_ref;Vel_o=Vel_i      # single phase neglect velocity change
        [DP_static,DP_hydrostatic,DP_dynamic,DP_total]=DP_ref.DP_ref_1ph(Re_ref,self.Ra_tube,rho_ref,G_ref_elem,Vel_i,Vel_o,z_i,z_o,self.L_elem,self.Dh_port)
        P_ri = P_ro-DP_total/1e3
        T_ri = CP.PropsSI('T','P',P_ri*1e3,'H',h_ri,self.fld)-273.15        # [C]
        h_ri_v=CP.PropsSI('H','P',P_ri*1e3,'Q',1,self.fld)            # [J/kg]
        h_ri_l=CP.PropsSI('H','P',P_ri*1e3,'Q',0,self.fld)
        x_ri = (h_ri-h_ri_l)/(h_ri_v-h_ri_l)
        
    else:           # TP
        T_ro = CP.PropsSI('T','P',P_ro*1e3,'H',h_ro,self.fld)-273.15        # [C]
        T_sat = T_ro
        
        q_H = 1.5e4         # [W/m^2]
        q_H_tmp = q_H+1e4
        while abs(q_H-q_H_tmp) > 10:
            q_H_tmp = q_H
            htc_ref = htc_2ph_Kim_2013(self.fld,q_H_tmp,self.G_ref_elem,T_sat,x_ro,self.Dh_port,self.P_tube_outer,self.P_port_tot,self.asp_ratio)
            htc_ref = htc_ref           # [W]
            while abs(Q_r-Q_a)>1e-3:
                h_f = CP.PropsSI('H','T',T_sa+273.15,'Q',0,'H2O')      # [J/kg]
                h_g = CP.PropsSI('H','T',T_sa+273.15,'Q',1,'H2O')
                h_fg = h_g - h_f
                if T_sa < Tdp_ai:
                    omega_s = CP.HAPropsSI('W','T',T_sa+273.15,'P',P_ai*1e3,'RH',1)             # [-]
                    Vda_s = CP.HAPropsSI('Vda','T',T_sa+273.15,'P',P_ai*1e3,'RH',1)             # [m^3/kg]
                    rho_As = omega_s / Vda_s            # [kg/m^3], kg of H2O/m^3
                    [Q_a_sen, Q_a_lat, Q_r, omega_Ao, T_ao] = dehumidification(T_ai,T_ro,x_ro,T_sa,htc_air,self.A_air_elem,cp_air,rho_air,Le_a,rho_As,\
                                                                        self.vol_dot_air_elem,Vda_i,m_dot_da,omega_Ai,h_fg,self.m_dot_ref_elem,cp_ref,self.R_wall,\
                                                                        htc_ref,self.A_ref_elem,rho_Ai,eta_t)
                    Q_a = Q_a_sen+Q_a_lat
                    Q = Q_a
                else:
                    omega_Ao = omega_Ai
                    C_ratio = 0
                    C_ref=999999999
                    UA = 1 / (R_air + 1/(htc_ref*self.A_ref_elem))
                    C_min = C_ha
                    NTU = UA/C_min
                    sigma = 1-np.exp(-NTU)
                    Q = sigma*C_min*(T_ai-T_ro)             # [kW]
                    Q_a = Q; Q_r = Q
                    Q_a_sen = Q_a; Q_a_lat = 0
                    T_ao = T_ai - Q/C_ha
                T_sa += (Q_a-Q_r)/10
            q_H = Q/self.A_ref_elem

        h_ri = h_ro - Q/m_dot_ref
        [DP_static,DP_hydrostatic,DP_dynamic,DP_total]=DP_ref.DP_ref_2ph(self,P_ro,self.G_ref_elem,x_ro,h_ri,z_i,z_o,self.L_elem,self.Dh_port,q_H)
        P_ri = P_ro-DP_total/1e3
        T_ri=CP.PropsSI('T','P',P_ri*1e3,'H',h_ri,self.fld)-273.15      # [C]
        h_ri_v=CP.PropsSI('H','P',P_ri*1e3,'Q',1,self.fld)      # [J/kg]
        h_ri_l=CP.PropsSI('H','P',P_ri*1e3,'Q',0,self.fld)      # [J/kg]
        x_ri = (h_ri-h_ri_l)/(h_ri_v-h_ri_l)
    
    SHR = Q_a_sen / Q_a*100         # [%], sensible heat ratio
    DELTAP_air_tmp = 0
    DELTAP_air = 20
    P_ao_tmp = self.P_ai - DELTAP_air_tmp
    while abs(DELTAP_air_tmp-DELTAP_air)>1e-8:
        DELTAP_air_tmp = DELTAP_air
        DELTAP_air = DP_air(self,self.T_ai,self.P_ai,T_ao,P_ao_tmp)
        P_ao = self.P_ai-DELTAP_air/1e3
    Tdp_ao = CP.HAPropsSI('Tdp','T',T_ao+273.15,'P',P_ao*1e3,'W',omega_Ao)-273.15           # [C]

    return Q,T_ri,T_ao,P_ri,x_ri,DP_total,htc_ref,h_ri,Q_a_sen,Q_a_lat,SHR,Tdp_ao,P_ao,DELTAP_air
import numpy as np
import CoolProp.CoolProp as CP
import simulation_utils.DP_ref as DP_ref
from simulation_utils.nusselt_1ph_gnielinskimodified import nusselt_1ph_gnielinskimodified
from simulation_utils.htc_2ph_chen_1966 import htc_2ph_chen_1966
from simulation_utils.htc_2ph_Kim_2013 import htc_2ph_Kim_2013

def epsilon_NTU(self,x_ro,h_ro,P_ro,m_dot_ref,T_ai,P_ai,m_dot_air,R_air,Re,G_ref_elem,htc_correlation,N_slab):
    '''
        epsilon-NTU method
        Since this is for cross flow and counter flow evaporator, we need to back iterate the whole process

         * Input Variables
            @x_ro           : refrigerant outlet vapor quality, [-]
            @h_ro           : refrigerant outlet enthalpy, [kJ/kg]
            @P_ro           : refrigerant outlet pressure, [kPa]
            @m_dot_ref      : refrigerant mass flow rate, [kg/s]
            @T_ai           : air inlet temperature, [C]
            @P_ai           : air inlet pressure, [kPa]
            @m_dot_air      : air mass flow rate, [kg/s]
            @R_air          : air thermal resistance, [K/W-m^2]
            @Re             : Reynolds number, [-]
            @G_ref_elem     : refrigerant mass velocity, [kg/m^2-s]
            @htc_correlation    : htc correlation method.
        
         * Output Variables
            @NTU, sigma
            @Q              : heat transfer rate, [kW]
            @T_ro           : refrigerant outlet temperature, [C]
            @T_ao           : air outlet temperature, [C]
    '''
    if N_slab%2 == 0:
        z_i = 0; z_o = self.L_elem
    else:
        z_i = self.L_elem; z_o = 0
    
    if (x_ro>1) or (x_ro<0):          # SP
        T_ro = CP.PropsSI('T','P',P_ro*1e3,'H',h_ro*1e3,self.fld)-273.15        # [C]
        Pr = CP.PropsSI('PRANDTL','T',T_ro+273.15,'P',P_ro*1e3,self.fld)        # [-]
        mu_ref = CP.PropsSI('V','T',T_ro+273.15,'P',P_ro*1e3,self.fld)              # [N/m^2-s]
        Re_ref = self.G_ref_elem*self.Dh_port/mu_ref
        Nu_1ph = nusselt_1ph_gnielinskimodified(Re_ref,Pr,self.asp_ratio)
        k_ref = CP.PropsSI('L','T',T_ro+273.15,'P',P_ro*1e3,self.fld)/1e3              # [kW/m-K]
        htc_ref = Nu_1ph*k_ref/self.Dh_port            # [kW/m^2-K]
        cp_ref = CP.PropsSI('C','T',T_ro+273.15,'P',P_ro*1e3,self.fld)/1e3        # [kJ/kg-K]
        rho_ref=CP.PropsSI('D','T',T_ro+273.15,'P',P_ro*1e3,self.fld)
        cp_air = CP.PropsSI('C','T',T_ai+273.15,'P',P_ai*1e3,'Air')/1e3
        C_ref = cp_ref*m_dot_ref            # [KJ/K]
        C_air = cp_air*m_dot_air
        C_min = min(C_ref,C_air)
        C_max = max(C_ref,C_air)
        C_ratio = C_min/C_max

        UA = 1 / (R_air + 1/(htc_ref*self.A_ref_elem))
        NTU = UA/C_min
        sigma = 1-np.exp((1/C_ratio)*NTU**0.22*(np.exp(-C_ratio*NTU**0.78)-1))

        T_ri = (C_ref*T_ro-sigma*C_min*T_ai)/(C_ref - sigma*C_min)
        Q = C_ref*(T_ro-T_ri)
        h_ri = h_ro - Q/m_dot_ref
        T_ao = T_ai - Q/C_air
        Vel_i=self.G_ref_elem/rho_ref;Vel_o=Vel_i      # single phase neglect velocity change
        [DP_static,DP_hydrostatic,DP_dynamic,DP_total]=DP_ref.DP_ref_1ph(Re,self.Ra_tube,rho_ref,G_ref_elem,Vel_i,Vel_o,z_i,z_o,self.L_elem,self.Dh_port)
        P_ri = P_ro-DP_total/1e3
        T_ri = CP.PropsSI('T','P',P_ri*1e3,'H',h_ri*1e3,self.fld)-273.15        # [C]
        h_ri_v=CP.PropsSI('H','P',P_ri*1e3,'Q',1,self.fld)/1e3            # [kJ/kg]
        h_ri_l=CP.PropsSI('H','P',P_ri*1e3,'Q',0,self.fld)/1e3 
        x_ri = (h_ri-h_ri_l)/(h_ri_v-h_ri_l)
        
    else:           # TP
        T_ro = CP.PropsSI('T','P',P_ro*1e3,'H',h_ro*1e3,self.fld)-273.15        # [C]
        T_sat = T_ro
        if htc_correlation == 'chen_1966':
            [htc_ref,T_wall_tmp] = htc_2ph_chen_1966(self.fld,T_ro,P_ro,x_ro,self.G_ref_elem,self.Dh_port,R_air,T_ai,self.A_ref_elem,self.asp_ratio)
            htc_ref = htc_ref/1e3
            C_ratio = 0
            C_ref=999999999
            cp_air = CP.PropsSI('C','T',T_ai+273.15,'P',P_ai*1e3,'Air')/1e3
            C_air = cp_air*m_dot_air
            UA = 1 / (R_air + 1/(htc_ref*self.A_ref_elem))
            C_min = C_air
            NTU = UA/C_min
            sigma = 1-np.exp(-NTU)
            Q = sigma*C_min*(T_ai-T_ro)
            q_H = Q/self.A_ref_elem*1e3
            
        if htc_correlation == 'kim_2013':
            q_H = 1.5e4         # [W/m^2]
            q_H_tmp = q_H+1e4
            while abs(q_H-q_H_tmp) > 10:
                q_H_tmp = q_H
                htc_ref = htc_2ph_Kim_2013(self.fld,q_H_tmp,self.G_ref_elem,T_sat,x_ro,self.Dh_port,self.P_tube_outer,self.P_port_tot,self.asp_ratio)
                htc_ref = htc_ref/1e3           # [W] -> [kW]
                C_ratio = 0
                C_ref=999999999
                cp_air = CP.PropsSI('C','T',T_ai+273.15,'P',P_ai*1e3,'Air')/1e3
                C_air = cp_air*m_dot_air
                UA = 1 / (R_air + 1/(htc_ref*self.A_ref_elem))
                C_min = C_air
                NTU = UA/C_min
                sigma = 1-np.exp(-NTU)

                Q = sigma*C_min*(T_ai-T_ro)             # [kW]
                q_H = Q/self.A_ref_elem*1e3
        T_ao = T_ai - Q/C_air
        h_ri = h_ro - Q/m_dot_ref
        [DP_static,DP_hydrostatic,DP_dynamic,DP_total]=DP_ref.DP_ref_2ph(self,P_ro,self.G_ref_elem,x_ro,h_ri,z_i,z_o,self.L_elem,self.Dh_port,q_H)
        P_ri = P_ro-DP_total/1e3
        T_ri=CP.PropsSI('T','P',P_ri*1e3,'H',h_ri*1e3,self.fld)-273.15      # [C]
        h_ri_v=CP.PropsSI('H','P',P_ri*1e3,'Q',1,self.fld)/1e3      # [kJ/kg]
        h_ri_l=CP.PropsSI('H','P',P_ri*1e3,'Q',0,self.fld)/1e3      # [kJ/kg]
        x_ri = (h_ri-h_ri_l)/(h_ri_v-h_ri_l)

    return NTU,sigma,Q,T_ri,T_ao,P_ri,x_ri,DP_total,htc_ref,h_ri
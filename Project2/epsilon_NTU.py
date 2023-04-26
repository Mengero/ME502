import numpy as np
import CoolProp.CoolProp as CP
import simulation_utils.DP_ref as DP_ref
from simulation_utils.htc_1ph import htc_1ph
from simulation_utils.htc_2ph_y1234yf_Shah import htc_2ph_y1234yf_shah
from simulation_utils.htc_2ph_r1234yf_Kim_2014 import htc_2ph_r1234yf_Kim_2014

def epsilon_NTU(self,x_ri,T_ri,P_ri,m_dot_ref,T_ai,P_ai,m_dot_air,R_air,Re,G_ref_elem,htc_correlation):
    '''
        epsilon-NTU method
        if the refrigerant is single phase, use equation (11.32)
        if the refrigerant is two phase, use equation (11.35a)

         * Input Variables
            @x_ref          : refrigerant inlet quality, [-]
            @htc_ref        : heat transfer coefficient of refrigerant, [kJ/m^2-K]
            @T_ri           : refrigerant inlet temperature, [C]
            @P_ri           : refrigerant inlet pressure, [kPa]
            @m_dot_ref      : refrigerant mass flow rate, [kg/s]
            @T_ai           : air inlet temperature, [C]
            @P_ai           : air inlet pressure, [kPa]
            @m_dot_air      : air mass flow rate, [kg/s]
            @R_air          : air side thermal resistance, [K/kJ]
        
         * Output Variables
            @NTU, sigma
            @Q              : heat transfer rate, [kW]
            @T_ro           : refrigerant outlet temperature, [C]
            @T_ao           : air outlet temperature, [C]
    '''
    
    if (x_ri>1) or (x_ri<0):          # SP
        h_ri = CP.PropsSI('H','T',T_ri+273.15,'P',P_ri*1e3,self.fld)/1e3        # [kJ/kg]
        htc_ref = htc_1ph(self,Re,T_ri,P_ri)/1e3            # [kW/m^2-K]
        cp_ref = CP.PropsSI('C','T',T_ri+273.15,'P',P_ri*1e3,self.fld)/1e3        # [kJ/kg-K]
        rho_ref=CP.PropsSI('D','T',T_ri+273.15,'P',P_ri*1e3,self.fld)
        cp_air = CP.PropsSI('C','T',T_ai+273.15,'P',P_ai*1e3,'Air')/1e3
        C_ref = cp_ref*m_dot_ref            # [KJ/K]
        C_air = cp_air*m_dot_air
        C_min = min(C_ref,C_air)
        C_max = max(C_ref,C_air)
        C_ratio = C_min/C_max

        UA = 1 / (R_air + 1/(htc_ref*self.A_ref_elem))
        NTU = UA/C_min
        sigma = 1-np.exp((1/C_ratio)*NTU**0.22*(np.exp(-C_ratio*NTU**0.78)-1))

        Q = sigma*C_min*(T_ri-T_ai)
        T_ao = T_ai + Q/C_air
        h_ro = h_ri - Q/m_dot_ref
        Vel_i=self.G_ref_elem/rho_ref;Vel_o=Vel_i      # single phase neglect velocity change
        z_i=0; z_o=0        # horizontal tube
        [DP_static,DP_hydrostatic,DP_dynamic,DP_total]=DP_ref.DP_ref_1ph(Re,self.Ra_tube,rho_ref,G_ref_elem,Vel_i,Vel_o,z_i,z_o,self.L_elem,self.Dh_port)
        P_ro = P_ri-DP_total/1e3
        T_ro = CP.PropsSI('T','P',P_ro*1e3,'H',h_ro*1e3,self.fld)-273.15        # [C]
        h_ro_v=CP.PropsSI('H','P',P_ro*1e3,'Q',1,self.fld)/1e3            # [kJ/kg]
        h_ro_l=CP.PropsSI('H','P',P_ro*1e3,'Q',0,self.fld)/1e3 
        x_ro = (h_ro-h_ro_l)/(h_ro_v-h_ro_l)
    else:           # TP
        h_ri = CP.PropsSI('H','P',P_ri*1e3,'Q',x_ri,self.fld)/1e3        # [kJ/kg]
        if htc_correlation == 'shah':
            htc_ref = htc_2ph_y1234yf_shah(self.fld,G_ref_elem,self.Dh_port,P_ri,x_ri)/1e3      # [kW/m^2-K]
        else:
            htc_ref = htc_2ph_r1234yf_Kim_2014(self.fld,G_ref_elem,self.Dh_port,self.asp_ratio,P_ri,x_ri)/1e3
        C_ratio = 0
        C_ref=999999999
        cp_air = CP.PropsSI('C','T',T_ai+273.15,'P',P_ai*1e3,'Air')/1e3
        C_air = cp_air*m_dot_air
        UA = 1 / (R_air + 1/(htc_ref*self.A_ref_elem))
        C_min = C_air
        NTU = UA/C_min
        sigma = 1-np.exp(-NTU)

        z_i=0;z_o=0
        Q = sigma*C_min*(T_ri-T_ai)
        T_ao = T_ai + Q/C_air
        h_ro = h_ri - Q/m_dot_ref
        [DP_static,DP_hydrostatic,DP_dynamic,DP_total]=DP_ref.DP_ref_2ph(self,Re,\
            self.Ra_tube,P_ri,G_ref_elem,x_ri,h_ro,z_i,z_o,self.L_elem,self.Dh_port)
        P_ro = P_ri-DP_total/1e3
        T_ro=CP.PropsSI('T','P',P_ro*1e3,'H',h_ro*1e3,self.fld)-273.15      # [C]
        h_ro_v=CP.PropsSI('H','P',P_ro*1e3,'Q',1,self.fld)/1e3      # [kJ/kg]
        h_ro_l=CP.PropsSI('H','P',P_ro*1e3,'Q',0,self.fld)/1e3      # [kJ/kg]
        x_ro = (h_ro-h_ro_l)/(h_ro_v-h_ro_l)

    return NTU,sigma,Q,T_ro,T_ao,P_ro,x_ro,DP_total,htc_ref
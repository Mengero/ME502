import numpy as np
import CoolProp.CoolProp as CP
from dehumidification_model import dehumidification
from simulation_utils.DP_air import DP_air
from simulation_utils.airsidehtc_chang_wang import airsidehtc_chang_wang
from simulation_utils.htc_H2OEG import htc_H2OEG
from simulation_utils.DP_H2OEG import DP_H2OEG

class cooling_coil:
    '''
        Unit system of the project: C, kJ, kPa, kW, s, kg ...
        Input Varibles:
            @x_ri               : [-], refrigerant inlet vapor quality
            @T_sat_ri           : [C], refrigerant inlet pressure
            @m_dot_ref          : [L/min], refrigerant volumetric flow rate
            @T_ai               : [C], air inlet temperature
            @P_ai               : [kPa], air inlet pressure
            @vol_dot_air        : [kg/h], air mass flow rate
            @N_segment          : [-], number of segment
    '''
    def __init__(self,T_ri,P_ri,vol_dot_r,T_ai,P_ai,vol_dot_a,RH_ai,N_segment):
        self.fld='INCOMP::MEG[0.5]'
        # -------------- Geometry Info ------------- #
        self.N_slab=1       # number of slabs
        self.N_pass=1       # number of passes
        self.L_tube_pass=300e-3             # [m], tube length in one pass
        self.N_tube=np.array([23])                 # [-], number of tubes in one pass
        self.N_tube_tot=np.sum(self.N_tube)
        self.t_wall=0.35e-3                 # [m], wall thickness
        self.t_tube=3e-3                    # [m], tube thickness, tube minor
        self.D_tube=35e-3                   # [m], tube depth, tube major
        self.N_port=11                       # [-], number of ports in one tube
        self.Ra_tube=1e-6                   # [m], roughness of tube inner surface
        self.theta_louver=15                # [deg], louver angle
        self.P_louver=1.3e-3                # [m], louver pitch
        self.L_louver=7.2e-3                # [m], louver length
        self.N_louverbank=2                 # [-], number of louver sets per fin
        self.h_fin=8e-3                     # [m], fin height
        self.t_fin=0.1e-3                   # [m], fin thickness
        self.P_fin=1.8e-3                   # [m], fin pitch (14 FPI)
        self.N_finptube = self.L_tube_pass/self.P_fin     # [-], number of fins on each tube
        self.D_fin=35e-3                    # [m], fin depth (same as tube)
        self.Dh_fin = 4*(self.P_fin-self.t_fin)*(self.h_fin-self.t_fin)/(((self.h_fin-self.t_fin)+(self.P_fin-self.t_fin))*2)

        self.H_port=self.t_tube-2*self.t_wall       # [m], height of port
        self.W_port=(self.D_tube-(self.N_port+1)*self.t_wall)/self.N_port       # [m], width of port
        self.P_port_tot = 2*(self.H_port+self.W_port)*(self.N_port-2)+2*(self.H_port+2*(self.W_port-self.H_port/2)+np.pi*self.H_port/2)       # [m], total perimeter of ports in one tube
        self.A_port_tot = (self.N_port-2)*(self.H_port*self.W_port)+2*((self.W_port-self.H_port/2)*self.H_port+np.pi*(self.H_port/2)**2/2)      # [m^2], total cross sectional area for refrigeratn in one tube
        self.Dh_port = 4*self.A_port_tot/self.P_port_tot
        self.asp_ratio = self.H_port/self.W_port
        self.P_tube_outer = (self.D_tube-self.t_tube)*2 + self.t_tube * np.pi     # [m], outer perimeter of one tube

        self.W_HX = self.L_tube_pass
        self.H_HX = self.t_tube*self.N_tube_tot+self.h_fin*(self.N_tube_tot-1)

        self.N_elemptube=N_segment
        self.N_elem = self.N_tube_tot*N_segment

        self.L_elem = self.L_tube_pass/self.N_elemptube
        self.P_tube = self.h_fin + self.t_tube
        
        # ------------- Area Calculation ----------- #
        self.A_ref_tot = self.P_port_tot*self.L_tube_pass*self.N_tube_tot
        self.A_ref_elem = self.A_ref_tot/self.N_elem
        self.A_fin_tot = (self.D_fin*self.h_fin)*2*self.N_finptube*self.N_tube_tot
        self.A_tube_tot = 2*(self.L_tube_pass-self.t_fin*self.N_finptube)*self.D_tube*(self.N_tube_tot-1)
        self.A_air_tot = self.A_fin_tot+self.A_tube_tot
        self.A_air_elem = self.A_air_tot/self.N_elem
        self.A_front = self.W_HX*self.H_HX      # frontal area of heat exchanger
        self.A_front_cross = self.L_tube_pass*self.t_tube*self.N_tube_tot + \
                       (self.t_fin*self.h_fin+self.t_fin*(self.P_fin-self.t_fin))*self.N_finptube*self.N_tube_tot # frontal cross sectional area of tube+fin
        self.A_air_free = self.A_front-self.A_front_cross           # minimum free flow area
        self.sigma = self.A_air_free/self.A_front
        self.k_w = 155/1e3      # [kW/m-K]
        self.rho_w = 2710   # [kg/m^3]
        self.R_wall = self.t_wall/(self.A_ref_elem*self.k_w)

        # ------------- Inlet Condition ------------ #
        self.T_ri = T_ri                                            # [C]
        self.P_ri = P_ri                                            # [kPa]
        self.vol_dot_r = vol_dot_r/60/1e3                           # [L/min] -> [m^3/s]
        self.vol_dot_a = vol_dot_a/2118.88                              # [CFM] -> [m^3/s]
        self.T_ai = T_ai                                                    # [C]
        self.P_ai = P_ai                                                    # [kPa]
        self.RH_ai = RH_ai                                                  # [-]
        self.rho_dai = 1/CP.HAPropsSI('Vda','T',self.T_ai+273.15,'P',self.P_ai*1e3,'RH',self.RH_ai)      # [kg/m^3]
        self.m_dot_da = self.vol_dot_a*self.rho_dai                       # [kg/s]
        self.m_dot_da_elem = self.m_dot_da/self.N_elem                    # [kg/s], dry air mass flow rate in each element
        self.rho_hai = 1/CP.HAPropsSI('Vha','T',self.T_ai+273.15,'P',self.P_ai*1e3,'RH',self.RH_ai)
        self.m_dot_ha = self.vol_dot_a*self.rho_hai
        self.m_dot_ha_elem = self.m_dot_ha/self.N_elem
        self.G_air = self.m_dot_da/self.A_air_free
        self.vel_air = self.vol_dot_a/self.A_air_free
        self.vol_dot_a_elem = self.vol_dot_a/self.N_elem
        self.rho_ri = CP.PropsSI('D','T',self.T_ri+273.15,'P',self.P_ri*1e3,self.fld)       # [kg/m^3]
        self.m_dot_ref = self.vol_dot_r*self.rho_ri
        self.m_dot_ref_elem = self.m_dot_ref/self.N_tube
        self.G_ref_elem = self.m_dot_ref_elem/self.A_port_tot
        self.vel_ref_elem = self.G_ref_elem/self.A_port_tot
        self.mu_ri = CP.PropsSI('V','T',self.T_ri+273.15,'P',self.P_ri*1e3,self.fld)        # [Pa-s]

def cc_tubeelem(self,T_ri,P_ri,T_ai,P_ai,RH_ai,R_air,m_dot_ref,m_dot_da,htc_air,eta_t,A_sa):
    '''
        Each tube element calculation
         * Input Variables:
            @T_ri           : refrigerant inlet temperature [C]
            @P_ri           : refrigerant inlet pressure [kPa]
        
         * Output Variables:
            @T_ro               : refrigerant outlet temperature [C]
            @P_ro               : refriegrant outlet pressure [kPa]
            @T_ao               : air outlet temperature [C]
            @P_ao               : air outlet pressure [kPa]
            @Q                  : heat transfer rate [kW]
            @DELTAP_air         : air side pressure drop [kPa]
            @DELTAP_H2OEG       : refrigerant side pressure drop [kPa]
    '''
    mu = CP.PropsSI('V','T',T_ri+273.15,'P',P_ri*1e3,self.fld)          # [Pa.s] 
    Re_ref = self.G_ref_elem * self.Dh_port / mu
    htc_ref = htc_H2OEG(Re_ref,T_ri,P_ri,self.W_port,self.H_port,self.Dh_port)          # [W/m^2-K]
    htc_ref = htc_ref*1e-3      # [kW/m^2-K]
    T_wall = (1 / R_air * T_ai + htc_ref * self.A_ref_elem * T_ri) / (htc_ref * self.A_ref_elem + 1 / R_air)
    h_ri = CP.PropsSI('H','T',T_ri+273.15,'P',P_ri*1e3,self.fld)/1e3        # [kJ/kg]
    rho_ref = CP.PropsSI('D','T',T_ri+273.15,'P',P_ri*1e3,self.fld)

    # Iteration with dehumidification model
    Q_r = 9999
    Q_a = 0
    Tdp_ai = CP.HAPropsSI('Tdp','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)-273.15
    T_sa = Tdp_ai-1

    cp_ref = CP.PropsSI('C','T',T_ri+273.15,'P',P_ri*1e3,self.fld)/1e3          # [kJ/kg-K]
    cp_air = CP.HAPropsSI('C','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)/1e3
    C_ref = cp_ref*m_dot_ref            # [KJ/K]
    C_air = cp_air*m_dot_da
    rho_air = 1/CP.HAPropsSI('Vha','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)          # [kg/m^3]
    k_air = CP.HAPropsSI('K','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)/1e3       # [kW/m-K]
    alpha_a = k_air / (rho_air*cp_air)               # [m^2/s]
    D_AB_a = 0.26e-4                                # [m^2/s]
    Le_a = alpha_a / D_AB_a
    Vda_i = CP.HAPropsSI('Vda','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)         # [m^3/kg]
    omega_Ai = CP.HAPropsSI('W','T',T_ai+273.15,'P',P_ai*1e3,'RH',RH_ai)        # [-]
    rho_Ai = omega_Ai / Vda_i

    # ------------------------- Energy Calculation (Included considerations for condensation) ----------------------- #
    while abs(Q_r-Q_a)>1e-7:
        h_f = CP.PropsSI('H','T',T_sa+273.15,'Q',0,'H2O')/1e3      # [kJ/kg]
        h_g = CP.PropsSI('H','T',T_sa+273.15,'Q',1,'H2O')/1e3
        h_fg = h_g - h_f
        if T_sa < Tdp_ai:
            omega_s = CP.HAPropsSI('W','T',T_sa+273.15,'P',P_ai*1e3,'RH',1)             # [-]
            Vda_s = CP.HAPropsSI('Vda','T',T_sa+273.15,'P',P_ai*1e3,'RH',1)             # [m^3/kg]
            rho_As = omega_s / Vda_s            # [kg/m^3], kg of H2O/m^3
            [Q_a_sen, Q_a_lat, Q_r, omega_Ao] = dehumidification(T_ai,T_ri,T_sa,htc_air,A_sa,cp_air,rho_air,Le_a,rho_As,\
                                    self.vol_dot_a_elem,Vda_i,self.m_dot_da_elem,omega_Ai,h_fg,self.m_dot_ref_elem,cp_ref,self.R_wall,\
                                    htc_ref,self.A_ref_elem,rho_Ai,eta_t)
            Q_a = Q_a_sen+Q_a_lat
            Q = Q_a
        else:
            C_min = min(C_ref,C_air)
            C_max = max(C_ref,C_air)
            C_ratio = C_min/C_max
            UA = 1 / (R_air + 1/(htc_ref*self.A_ref_elem))
            NTU = UA/C_min
            sigma = 1-np.exp((1/C_ratio)*NTU**0.22*(np.exp(-C_ratio*NTU**0.78)-1))

            Q = sigma*C_min*(T_ai-T_ri)
            Q_r = Q; Q_a = Q
            Q_a_sen = Q; Q_a_lat = 0

        if Q_a>Q_r:
            T_sa += (Q_a-Q_r)*10
        elif Q_a<Q_r:
            T_sa += (Q_a-Q_r)*10
        
    # ------------------------- Output Calculation ----------------------- #
    T_ao = T_ai - Q/C_air
    DELTAP_air_tmp = 0
    DELTAP_air = 20
    P_ao_tmp = P_ai - DELTAP_air_tmp/1e3
    while abs(DELTAP_air_tmp-DELTAP_air)>1e-5:
        DELTAP_air_tmp = DELTAP_air
        DELTAP_air = DP_air(self.sigma,T_ai,P_ai,T_ao,P_ao_tmp,self.vel_air,self.P_louver,self.P_fin,self.h_fin,self.t_fin,self.Dh_fin,self.L_louver,self.D_fin,\
           self.P_tube,self.t_tube,self.theta_louver,self.A_air_free,self.A_air_tot,self.A_front,self.G_air)             # [Pa]
        P_ao = P_ai-DELTAP_air/1e3                                              # [kPa]
    Tdp_ao = CP.HAPropsSI('Tdp','T',T_ao+273.15,'P',P_ao*1e3,'W',omega_Ao)-273.15           # [C]
    T_ro = T_ri + Q/C_ref
    DP_r = DP_H2OEG(Re_ref,self.Ra_tube,rho_ref,self.G_ref_elem,self.vel_ref_elem,self.vel_ref_elem,0,0,self.L_elem,self.Dh_port)
    P_ro = P_ri-DP_r/1e3            # [kPa]

    return T_ro,T_ai,Tdp_ai,T_ao,Tdp_ao,P_ro,htc_ref,Q_r,Q_a_sen,Q_a_lat

def cc_tube(self,T_ri,P_ri,T_ai,P_ai,RH_ai):
    variables = ['EG Temperature','air inlet temperature','air inlet dew point temperature',\
                  'air outlet temperature','air outlet dew point temperature','EG Pressure',\
                  'EG HTC','segment heat transfer rate','segment sensible heat transfer rate',\
                  'segment latent heat transfer rate']
    data = {variable: [] for variable in variables}
    x_data = range(1,41)

    [htc_air, eta_f, eta_t] = airsidehtc_chang_wang(P_ai,self.k_w,self.theta_louver,self.P_fin,self.P_louver,\
        self.h_fin,self.D_fin,self.L_louver,self.D_tube,self.t_tube,self.t_fin,T_ai,self.vel_air)
    htc_air = htc_air/1e3
    R_air = 1 / (htc_air * self.A_air_elem * eta_t)

    for i in range(self.N_elemptube):
        print(i)
        seg_result_tmp = cc_tubeelem(self,T_ri,P_ri,T_ai,P_ai,RH_ai,R_air,self.m_dot_ref_elem,self.m_dot_da_elem,htc_air,eta_t,self.A_air_elem)
        for i, value in enumerate(seg_result_tmp):
            data[variables[i]].append(value)
        T_ri = seg_result_tmp[0]
        P_ri = seg_result_tmp[5]
        
    return data,variables


def cc_system(self,T_ri,P_ri,T_ai,P_ai,RH_ai):
    '''
        Calculate the condenser system
        single slab, each tube identical
    '''
    [data,variables] = cc_tube(self,T_ri,P_ri,T_ai,P_ai,RH_ai)
    T_ro = data[variables[0]][-1]
    T_ao = np.mean(data[variables[3]])
    Tdp_ao = np.mean(data[variables[4]])
    P_ro = data[variables[5]][-1]
    Q_seg_tot = self.N_tube_tot * sum(data[variables[7]])
    Q_sen_tot = self.N_tube_tot * sum(data[variables[8]])
    Q_lat_tot = self.N_tube_tot * sum(data[variables[9]])

    return  data,variables,T_ro,P_ro,T_ao,Tdp_ao,Q_seg_tot,Q_sen_tot,Q_lat_tot
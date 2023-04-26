import numpy as np
import CoolProp.CoolProp as CP
from simulation_utils.DP_air import DP_air
from simulation_utils.airsidehtc_chang_wang import airsidehtc_chang_wang
from epsilon_NTU import epsilon_NTU

class condenser:
    '''
        Unit system of the project: C, kJ, kPa, kW, s, kg ...
        Input Varibles:
            @T_ri               : [C], refrigerant inlet temperature
            @T_sat_ri           : [C], refrigerant inlet pressure
            @m_dot_ref          : [L/min], refrigerant volumetric flow rate
            @T_ai               : [C], air inlet temperature
            @P_ai               : [kPa], air inlet pressure
            @vol_dot_air        : [kg/h], air mass flow rate
            @N_segment          : [-], number of segment
    '''
    def __init__(self,T_ri,T_sat_ri,m_dot_ref,T_ai,P_ai,vol_dot_air,N_segment):
        self.fld='R1234yf'
        # -------------- Geometry Info ------------- #
        self.N_slab=1       # number of slabs
        self.N_pass=4       # number of passes
        self.L_tube_pass=650e-3             # [m], tube length in one pass
        self.N_tube=np.array([20,15,10,5])                 # [-], number of tubes in one pass
        self.N_tube_tot=np.sum(self.N_tube)
        self.t_wall=0.35e-3                 # [m], wall thickness
        self.t_tube=1.7e-3                    # [m], tube thickness, tube minor
        self.D_tube=20e-3                   # [m], tube depth, tube major
        self.N_port=13                       # [-], number of ports in one tube
        self.Ra_tube=1e-6                   # [m], roughness of tube inner surface
        self.theta_louver=20                # [deg], louver angle
        self.P_louver=1.3e-3                # [m], louver pitch
        self.L_louver=7.2e-3                # [m], louver length
        self.N_louverbank=2                 # [-], number of louver sets per fin
        self.h_fin=8e-3                     # [m], fin height
        self.t_fin=0.1e-3                   # [m], fin thickness
        self.P_fin=1.4e-3                   # [m], fin pitch (18 FPI)
        self.N_finptube = self.L_tube_pass/self.P_fin     # [-], number of fins on each tube
        self.D_fin=20e-3                    # [m], fin depth (same as tube)
        self.Dh_fin = 4*(self.P_fin-self.t_fin)*(self.h_fin-self.t_fin)/(((self.h_fin-self.t_fin)+(self.P_fin-self.t_fin))*2)

        self.H_port=self.t_tube-2*self.t_wall       # [m], height of port
        self.W_port=(self.D_tube-(self.N_port+1)*self.t_wall)/self.N_port       # [m], width of port
        self.P_port_tot = 2*(self.H_port+self.W_port)*(self.N_port-2)+2*(self.H_port+2*(self.W_port-self.H_port/2)+np.pi*self.H_port/2)       # [m], total perimeter of ports in one tube
        self.A_port_tot = (self.N_port-2)*(self.H_port*self.W_port)+2*((self.W_port-self.H_port/2)*self.H_port+np.pi*(self.H_port/2)**2/2)      # [m^2], total cross sectional area for refrigeratn in one tube
        self.Dh_port = 4*self.A_port_tot/self.P_port_tot
        self.asp_ratio = self.H_port/self.W_port

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
        self.k_w = 155      # [W/m-K]
        self.rho_w = 2710   # [kg/m^3]

        # ------------- Inlet Condition ------------ #
        self.T_ri = T_ri                                                    # [C]
        self.T_sat_ri = T_sat_ri                                            # [C]
        self.m_dot_ref = m_dot_ref                                          # [kg/s]
        self.vol_dot_air = vol_dot_air/2118.88                              # [m^3/s]
        self.T_ai = T_ai                                                    # [C]
        self.P_ai = P_ai                                                    # [kPa]
        self.rho_ai = CP.PropsSI('D','T',self.T_ai+273.15,'P',self.P_ai*1e3,'Air')      # [kg/m^3]
        self.m_dot_air = self.vol_dot_air*self.rho_ai                       # [kg/s]
        self.m_dot_air_elem = self.m_dot_air/self.N_elem                    # [kg/s]
        self.G_air = self.m_dot_air/self.A_air_free
        self.vel_air = self.m_dot_air/self.rho_ai/self.A_air_free
        self.vel_air = 4.758
        self.P_ri = CP.PropsSI('P','T',T_sat_ri+273.15,'Q',1,self.fld)/1e3      # [kPa]
        self.h_ri = CP.PropsSI('H','P',self.P_ri*1e3,'T',T_ri+273.15,self.fld)/1e3  # [kJ/kg]
        self.h_ri_f = CP.PropsSI('H','P',self.P_ri*1e3,'Q',0,self.fld)/1e3  # [kJ/kg]
        self.h_ri_g = CP.PropsSI('H','P',self.P_ri*1e3,'Q',1,self.fld)/1e3  # [kJ/kg]
        self.x_ri = (self.h_ri-self.h_ri_f)/(self.h_ri_g-self.h_ri_f)
        self.m_dot_ref_elem = self.m_dot_ref/self.N_tube
        self.G_ref_elem = self.m_dot_ref_elem/self.A_port_tot

def condenser_tubeelem(self,T_ri,P_ri,x_ri,R_air,index_pass,htc_correlation):
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
    if x_ri<0 or x_ri>1:            # SP
        rho_ri = CP.PropsSI('D','T',T_ri+273.15,'P',P_ri*1e3,self.fld)
        mu_ri = CP.PropsSI('V','T',T_ri+273.15,'P',P_ri*1e3,self.fld)
        Re_ref = self.G_ref_elem[index_pass]*self.Dh_port/mu_ri
        [NTU,sigma,Q,T_ro,T_ao,P_ro,x_ro,DELTAP_ref,htc_ref] = epsilon_NTU(self,x_ri,T_ri,P_ri,\
            self.m_dot_ref_elem[index_pass],self.T_ai,self.P_ai,self.m_dot_air_elem,R_air,Re_ref,self.G_ref_elem[index_pass],htc_correlation)
    else:                           # TP
        rho_ri = CP.PropsSI('D','P',P_ri*1e3,'Q',x_ri,self.fld)
        mu_ri = CP.PropsSI('V','P',P_ri*1e3,'Q',x_ri,self.fld)
        Re_ref = self.G_ref_elem[index_pass]*self.Dh_port/mu_ri
        [NTU,sigma,Q,T_ro,T_ao,P_ro,x_ro,DELTAP_ref,htc_ref] = epsilon_NTU(self,x_ri,T_ri,P_ri,\
            self.m_dot_ref_elem[index_pass],self.T_ai,self.P_ai,self.m_dot_air_elem,R_air,Re_ref,self.G_ref_elem[index_pass],htc_correlation)

    DELTAP_air_tmp = 0
    DELTAP_air = 20
    P_ao_tmp = self.P_ai - DELTAP_air_tmp
    while abs(DELTAP_air_tmp-DELTAP_air)>1e-8:
        DELTAP_air_tmp = DELTAP_air
        DELTAP_air = DP_air(self,self.T_ai,self.P_ai,T_ao,P_ao_tmp)
        P_ao = self.P_ai-DELTAP_air/1e3

    return T_ro,P_ro,x_ro,T_ao,P_ao,Q,DELTAP_air,DELTAP_ref,htc_ref

def condenser_tube(self,T_ri,P_ri,x_ri,R_air,index_pass,htc_correlation):
    T_r_dist = np.array([])
    P_r_dist = np.array([])
    x_r_dist = np.array([])
    T_a_dist = np.array([])
    P_a_dist = np.array([])
    Q_dist = np.array([])
    DELTAP_air_dist = np.array([])
    DELTAP_ref_dist = np.array([])
    htc_ref_dist = np.array([])

    if self.N_elemptube==1:         # Lumped Condition
        x_r = 0.5
        [T_ro,P_ro,x_ro,T_ao,P_ao,Q,DELTAP_air,DELTAP_ref,htc_ref]=condenser_tubeelem(self,T_ri,P_ri,x_ri,R_air,index_pass,htc_correlation)
        T_r_dist = np.append(T_r_dist,T_ro)
        T_a_dist = np.append(T_a_dist,T_ao)
        P_r_dist = np.append(P_r_dist,P_ro)
        P_a_dist = np.append(P_a_dist,P_ao)
        x_r_dist = np.append(x_r_dist,x_ro)
        Q_dist = np.append(Q_dist,Q)
        DELTAP_air_dist = np.append(DELTAP_air_dist,DELTAP_air)
        DELTAP_ref_dist = np.append(DELTAP_ref_dist,DELTAP_ref)
        htc_ref_dist = np.append(htc_ref_dist,htc_ref)
        T_ri = T_ro
        P_ri = P_ro
        x_ri = x_ro

    else:
        for i in range(self.N_elemptube):
            [T_ro,P_ro,x_ro,T_ao,P_ao,Q,DELTAP_air,DELTAP_ref,htc_ref]=condenser_tubeelem(self,T_ri,P_ri,x_ri,R_air,index_pass,htc_correlation)
            T_r_dist = np.append(T_r_dist,T_ro)
            T_a_dist = np.append(T_a_dist,T_ao)
            P_r_dist = np.append(P_r_dist,P_ro)
            P_a_dist = np.append(P_a_dist,P_ao)
            x_r_dist = np.append(x_r_dist,x_ro)
            Q_dist = np.append(Q_dist,Q)
            DELTAP_air_dist = np.append(DELTAP_air_dist,DELTAP_air)
            DELTAP_ref_dist = np.append(DELTAP_ref_dist,DELTAP_ref)
            htc_ref_dist = np.append(htc_ref_dist,htc_ref)
            T_ri = T_ro
            P_ri = P_ro
            x_ri = x_ro
        
    return T_ro,P_ro,x_ro,T_r_dist,x_r_dist,T_a_dist,P_r_dist,P_a_dist,Q_dist,DELTAP_air_dist,DELTAP_ref_dist,htc_ref_dist


def condenser_system(self,htc_correlation):
    '''
        Calculate the condenser system
    '''
    # First determine the heat transfer coefficient of air (same for all the elements)
    [htc_air, eta_f, eta_t] = airsidehtc_chang_wang(self.P_ai,self.k_w,self.theta_louver,self.P_fin,self.P_louver,\
        self.h_fin,self.D_fin,self.L_louver,self.D_tube,self.t_tube,self.t_fin,self.T_ai,self.vel_air)
    self.htc_air = htc_air/1e3
    self.R_air = 1 / (self.htc_air * self.A_air_elem * eta_t)
    T_ro_dict = {}
    P_ro_dict = {}
    x_ro_dict = {}
    x_r_dist_dict = {}
    T_r_dist_dict = {}
    T_a_dist_dict = {}
    P_r_dist_dict = {}
    P_a_dist_dict = {}
    Q_dist_dict = {}
    DELTAP_air_dist_dict = {}
    DELTAP_ref_dist_dict = {}
    htc_ref_dist_dict = {}

    P_ri = self.P_ri
    T_ri = self.T_ri
    x_ri = self.x_ri
    # Since uniform distribution, only calculate one tube and leads to all the other tubes
    # Iteration twice for each pass seperatedly
    for index_pass in range(self.N_pass):
        [T_ro,P_ro,x_ro,T_r_dist,x_r_dist,T_a_dist,P_r_dist,P_a_dist,Q_dist,DELTAP_air_dist,DELTAP_ref_dist,htc_ref_dist] = condenser_tube(self,T_ri,P_ri,x_ri,self.R_air,index_pass,htc_correlation)
        T_ro_dict[index_pass] = T_ro
        P_ro_dict[index_pass] = P_ro
        x_ro_dict[index_pass] = x_ro
        T_r_dist_dict[index_pass] = T_r_dist
        T_a_dist_dict[index_pass] = T_a_dist
        P_r_dist_dict[index_pass] = P_r_dist
        P_a_dist_dict[index_pass] = P_a_dist
        x_r_dist_dict[index_pass] = x_r_dist
        Q_dist_dict[index_pass] = Q_dist
        DELTAP_air_dist_dict[index_pass] = DELTAP_air_dist
        DELTAP_ref_dist_dict[index_pass] = DELTAP_ref_dist
        htc_ref_dist_dict[index_pass] = htc_ref_dist
        T_ri = T_ro
        P_ri = P_ro
        x_ri = x_ro

    T_ao = (np.average(T_a_dist_dict[0]) + np.average(T_a_dist_dict[1]))/2
    P_ao = (np.average(P_a_dist_dict[0]) + np.average(P_a_dist_dict[1]))/2

    Q_tot = np.sum(Q_dist_dict[0])*self.N_tube[0] + np.sum(Q_dist_dict[1])*self.N_tube[1]

    return  T_ro,P_ro,x_ro,T_ao,P_ao,Q_tot,T_ro_dict,P_ro_dict,x_ro_dict,Q_dist_dict,T_r_dist_dict,T_a_dist_dict,P_r_dist_dict,P_a_dist_dict,DELTAP_air_dist_dict,DELTAP_ref_dist_dict,x_r_dist_dict,htc_ref_dist_dict


# -------------------- Test Session ---------------------- #

# fld = 'R1234yf'
# T_ri = 75       # [C]
# T_sat_ri = 48   # [C]
# m_dot_cr = 35/1e3   # [g/s]
# T_ai = 35       # [C]
# P_ai = 99.5     # [kPa]
# vol_dot_cai = 1500    # [CFM], cubic feet per minute
# N_segment = 60
# self=condenser(T_ri,T_sat_ri,m_dot_cr,T_ai,P_ai,vol_dot_cai,N_segment)
# htc_correlation = 'Kim'

# [T_ro,P_ro,x_ro,T_ao,P_ao,Q_tot,T_ro_dict,P_ro_dict,x_ro_dict,Q_dist_dict,T_r_dist_dict,T_a_dist_dict,P_r_dist_dict,P_a_dist_dict,DELTAP_air_dist_dict,DELTAP_ref_dist_dict,x_r_dist_dict,htc_ref_dist_dict] = condenser_system(self,htc_correlation)
# print(T_ro)
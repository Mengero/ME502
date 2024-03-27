import numpy as np
import CoolProp.CoolProp as CP
from simulation_utils.DP_air import DP_air
from simulation_utils.airsidehtc_chang_wang import airsidehtc_chang_wang
from epsilon_NTU import epsilon_NTU
from dehumidification_model import dehumidification
from simulation_utils.nusselt_1ph_gnielinskimodified import nusselt_1ph_gnielinskimodified
from simulation_utils.htc_2ph_Kim_2013 import htc_2ph_Kim_2013

class evaporator:
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
    def __init__(self,h_ri,Tsh_ro,m_dot_ref,T_ai,P_ai,RH_ai,m_dot_air,N_segment):
        self.fld='R1234yf'
        # -------------- Geometry Info ------------- #
        self.N_slab=4       # number of slabs
        self.N_pass=1       # number of passes
        self.L_tube_pass=300e-3             # [m], tube length in one pass
        self.N_tube=np.array([25])                 # [-], number of tubes in one pass
        self.N_tube_tot=np.sum(self.N_tube)
        self.t_wall=0.35e-3                 # [m], wall thickness
        self.t_tube=1.7e-3                    # [m], tube thickness, tube minor
        self.D_tube=10e-3                   # [m], tube depth, tube major
        self.N_port=7                       # [-], number of ports in one tube
        self.Ra_tube=1e-6                   # [m], roughness of tube inner surface
        self.theta_louver=15                # [deg], louver angle
        self.P_louver=1.3e-3                # [m], louver pitch
        self.L_louver=7.2e-3                # [m], louver length
        self.N_louverbank=2                 # [-], number of louver sets per fin
        self.h_fin=8e-3                     # [m], fin height
        self.t_fin=0.1e-3                   # [m], fin thickness
        self.P_fin=1.8e-3                   # [m], fin pitch (14 FPI)
        self.N_finptube = self.L_tube_pass/self.P_fin     # [-], number of fins on each tube
        self.D_fin=10e-3                    # [m], fin depth (same as tube)
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
        self.A_air_free_elem = self.A_air_free/self.N_elem
        self.sigma = self.A_air_free/self.A_front
        self.k_w = 155      # [W/m-K]
        self.rho_w = 2710   # [kg/m^3]

        # ------------- Inlet Condition ------------ #
        self.h_ri = h_ri                                                    # [kJ/kg]
        self.Tsh_ro = Tsh_ro                                                # [Cs]
        self.m_dot_ref = m_dot_ref/1e3                                      # [g/s] -> [kg/s]
        self.m_dot_air = m_dot_air/60                                       # [kg/min] -> [kg/s]
        self.T_ai = T_ai                                                    # [C]
        self.P_ai = P_ai                                                    # [kPa]
        self.RH_ai = RH_ai                                                  # [-]
        self.rho_ai = 1/CP.HAPropsSI('V','T',self.T_ai+273.15,'P',self.P_ai*1e3,'RH',RH_ai)     # [kg/m^3]
        self.vol_dot_air = self.m_dot_air/self.rho_ai                            # [m^3/a]
        self.vol_dot_air_elem = self.vol_dot_air/self.N_elem
        self.m_dot_air_elem = self.m_dot_air/self.N_elem                    # [kg/s]
        self.G_air = self.m_dot_air/self.A_air_free
        self.vel_air = self.m_dot_air/self.rho_ai/self.A_air_free
        self.m_dot_ref_elem = self.m_dot_ref/self.N_tube_tot
        self.G_ref_elem = self.m_dot_ref_elem/self.A_port_tot
        self.Tdp_ai = CP.HAPropsSI('Tdp','T',self.T_ai+273.15,'P',self.P_ai*1e3,'RH',RH_ai)-273.15

        self.R_wall = self.t_wall/(self.A_ref_elem*self.k_w)

def evaporator_tubeelem(self,h_ro,P_ro,x_ro,T_ai,P_ai,RH_ai,N_slab):
    '''
        Each tube element calculation
         * Input Variables:
            @h_ro               : refrigerant outlet enthalpy, [J/kg]
            @P_ro               : refrigerant outlet pressure, [kPa]
            @x_ro               : refrigerant outlet quality, [-]
            @T_ai               : air inlet tempreature, [C]
            @P_ai               : air inlet pressure, [kPa]
            @RH_ai              : air inlet relative humidity, [-]
        
         * Output Variables:
            @T_ro               : refrigerant outlet temperature [C]
            @P_ro               : refriegrant outlet pressure [kPa]
            @T_ao               : air outlet temperature [C]
            @P_ao               : air outlet pressure [kPa]
            @Q                  : heat transfer rate [kW]
            @DELTAP_air         : air side pressure drop [kPa]
            @DELTAP_H2OEG       : refrigerant side pressure drop [kPa]
    '''
    # ---------------- Air property calculation ----------------- #
    rho_ai = 1/CP.HAPropsSI('V','T',self.T_ai+273.15,'P',self.P_ai*1e3,'RH',RH_ai)     # [kg/m^3]
    m_dot_air = self.vol_dot_air_elem*rho_ai
    vel_air = self.vol_dot_air_elem/self.A_air_free_elem                # [m/s]
    rho_dai = 1/CP.HAPropsSI('Vda','T',self.T_ai+273.15,'P',self.P_ai*1e3,'RH',RH_ai)     # [kg/m^3]
    m_dot_da = self.vol_dot_air_elem*rho_dai
    [htc_air, eta_f, eta_t] = airsidehtc_chang_wang(P_ai,self.k_w,self.theta_louver,self.P_fin,self.P_louver,\
        self.h_fin,self.D_fin,self.L_louver,self.D_tube,self.t_tube,self.t_fin,T_ai,vel_air)
    htc_air = htc_air
    R_air = 1 / (htc_air * self.A_air_elem * eta_t)
    [Q,T_ri,T_ao,P_ri,x_ri,DELTAP_ref,htc_ref,h_ri,Q_a_sen,Q_a_lat,SHR,Tdp_ao,P_ao,DELTAP_air] = \
        epsilon_NTU(self,x_ro,h_ro,P_ro,self.m_dot_ref_elem,T_ai,P_ai,RH_ai,m_dot_air,m_dot_da,R_air,htc_air,eta_t,self.G_ref_elem,N_slab)

    return T_ri,P_ri,x_ri,h_ri,T_ao,P_ao,Tdp_ao,Q,DELTAP_air,DELTAP_ref,htc_ref,Q_a_sen,Q_a_lat,SHR

def evaporator_tube(self,T_ro,h_ro,P_ro,x_ro,T_air_dist,P_air_dist,Tdp_air_dist,N_slab):
    T_r_dist = np.zeros(self.N_elemptube)
    P_r_dist = np.zeros(self.N_elemptube)
    x_r_dist = np.zeros(self.N_elemptube)
    T_a_dist = np.zeros(self.N_elemptube)
    P_a_dist = np.zeros(self.N_elemptube)
    Q_dist = np.zeros(self.N_elemptube)
    SHR_dist = np.zeros(self.N_elemptube)
    Tdp_a_dist = np.zeros(self.N_elemptube)
    Q_a_sen_dist = np.zeros(self.N_elemptube)
    Q_a_lat_dist = np.zeros(self.N_elemptube)
    DELTAP_air_dist = np.zeros(self.N_elemptube)
    DELTAP_ref_dist = np.zeros(self.N_elemptube)
    htc_ref_dist = np.zeros(self.N_elemptube)

    RH_air_dist = CP.HAPropsSI('RH','T',T_air_dist+273.15,'P',P_air_dist*1e3,'Tdp',Tdp_air_dist+273.15)

    for i in range(self.N_elemptube):
        i = self.N_elemptube-1-i            # j follows the coordinate, i follows the flow direction of refrigerant
        [T_ri,P_ri,x_ri,h_ri,T_ao,P_ao,Tdp_ao,Q,DELTAP_air,DELTAP_ref,htc_ref,Q_a_sen,Q_a_lat,SHR]=evaporator_tubeelem(self,h_ro,P_ro,x_ro,T_air_dist[i],P_air_dist[i],RH_air_dist[i],N_slab)
        T_r_dist[i] = T_ri
        T_a_dist[i] = T_ao
        P_r_dist[i] = P_ri
        P_a_dist[i] = P_ao
        x_r_dist[i] = x_ri
        Q_dist[i] = Q
        Q_a_sen_dist[i] = Q_a_sen
        Q_a_lat_dist[i] = Q_a_lat
        SHR_dist[i] = SHR
        Tdp_a_dist[i] = Tdp_ao
        DELTAP_air_dist[i] = DELTAP_air
        DELTAP_ref_dist[i] = DELTAP_ref
        htc_ref_dist[i] = htc_ref
        T_ro = T_ri
        P_ro = P_ri
        x_ro = x_ri
        h_ro = h_ri
        
    return T_ri,P_ri,x_ri,h_ri,T_r_dist,x_r_dist,T_a_dist,P_r_dist,P_a_dist,Q_dist,DELTAP_air_dist,DELTAP_ref_dist,htc_ref_dist,Q_a_sen_dist,Q_a_lat_dist,SHR_dist,Tdp_a_dist

def result_organizer(result_tmp):
    '''
    variables in result_tmp:
        T_ri,P_ri,x_ri,T_ro,P_ro,x_ro,T_ao,P_ao,Q_tot,T_ro_dict,P_ro_dict,
        x_ro_dict,Q_dist_dict,Q_a_sen_dist_dict,Q_a_lat_dist_dict,T_r_dist_dict,
        T_a_dist_dict,P_r_dist_dict,P_a_dist_dict,DELTAP_air_dist_dict,
        DELTAP_ref_dist_dict,x_r_dist_dict,htc_ref_dist_dict,Tdp_ao_dist_dict,SHR_dist,T_ai_dist
    '''
    result = {}
    variables = ['T_ri','P_ri','x_ri','T_ro','P_ro','x_ro','T_ao','P_ao','Q_tot','T_ro_dist',\
                 'P_ro_dist','x_ro_dist','Q_dist','Q_a_sen_dist','Q_a_lat_dist','T_r_dist','T_a_dist',\
                 'P_r_dist','P_a_dist','DP_a_dist','DP_r_dist','x_r_dist','htc_r_dist','Tdp_a_dist',\
                 'SHR_dist','T_ai_dist']
    for i, value in enumerate(result_tmp):
        if i<=11:
            result[variables[i]] = (value)
        else:
            result[variables[i]] = np.concatenate((value[1],value[2],value[3],value[4]))

    return result

def evaporator_system(Ref,T_cro,P_cro,m_dot_ref,Tsh_ro,T_ai,P_ai,RH_ai,m_dot_air,N_segment):
    '''
        Calculate the condenser system
    '''
    h_ri = CP.PropsSI('H','T',T_cro+273.15,'P',P_cro*1e3,Ref)
    self = evaporator(h_ri,Tsh_ro,m_dot_ref,T_ai,P_ai,RH_ai,m_dot_air,N_segment)

    T_ro_dict = {}
    P_ro_dict = {}
    x_ro_dict = {}
    x_r_dist_dict = {}
    T_r_dist_dict = {}
    T_a_dist_dict = {}
    P_r_dist_dict = {}
    P_a_dist_dict = {}
    Q_dist_dict = {}
    SHR_dist_dict = {}
    Tdp_ao_dist_dict = {}
    Q_a_sen_dist_dict = {}
    Q_a_lat_dist_dict = {}
    DELTAP_air_dist_dict = {}
    DELTAP_ref_dist_dict = {}
    htc_ref_dist_dict = {}
    # Since uniform distribution, only calculate one tube and leads to all the other tubes
    # Iteration twice for each pass seperatedly
    P_ro_guess = P_cro*0.6
    P_ro_converge = np.array([])
    h_ri_tmp = h_ri + 100            # guessed value
    iteration_count = 0
    while abs(h_ri_tmp - h_ri) > 50:
        iteration_count+=1
        print('itertaions count for {}'.format(iteration_count))
        print(abs(h_ri-h_ri_tmp))
        P_ro_tmp = P_ro_guess               # regive determined guessed output pressure
        P_ro_converge = np.append(P_ro_converge, P_ro_tmp)
        T_ro_sat = CP.PropsSI('T','P',P_ro_tmp*1e3,'Q',0.5,self.fld)-273.15
        T_ro_tmp = Tsh_ro + T_ro_sat
        h_ro_tmp = CP.PropsSI('H','P',P_ro_tmp*1e3,'T',T_ro_tmp+273.15,self.fld)         # [J/kg]
        h_ro_v_tmp = CP.PropsSI('H','P',P_ro_tmp*1e3,'Q',1,self.fld)        # [J/kg]
        h_ro_l_tmp = CP.PropsSI('H','P',P_ro_tmp*1e3,'Q',0,self.fld)
        x_ro_tmp = (h_ro_tmp-h_ro_l_tmp)/(h_ro_v_tmp-h_ro_l_tmp)
        # Define air side condition every time before iteration
        T_ai_dist = {}
        P_ai_dist = {}
        Tdp_ai_dist = {}
        T_ai_dist[self.N_slab] = np.ones(self.N_elemptube)*self.T_ai
        P_ai_dist[self.N_slab] = np.ones(self.N_elemptube)*self.P_ai
        Tdp_ai_dist[self.N_slab] = np.ones(self.N_elemptube)*self.Tdp_ai
        # Iteration through each slab
        for N_slab in range(self.N_slab,0,-1):
            [T_ri_tmp,P_ri_tmp,x_ri_tmp,h_ri_tmp,T_r_dist,x_r_dist,T_a_dist,P_r_dist,P_a_dist,Q_dist,DELTAP_air_dist,DELTAP_ref_dist,\
             htc_ref_dist,Q_a_sen_dist,Q_a_lat_dist,SHR_dist,Tdp_a_dist] = evaporator_tube(self,T_ro_tmp,h_ro_tmp,P_ro_tmp,x_ro_tmp,\
                                                                                            T_ai_dist[N_slab],P_ai_dist[N_slab],Tdp_ai_dist[N_slab],N_slab)
            T_ro_dict[N_slab] = T_ro_tmp
            P_ro_dict[N_slab] = P_ro_tmp
            x_ro_dict[N_slab] = x_ro_tmp
            T_r_dist_dict[N_slab] = T_r_dist
            T_a_dist_dict[N_slab] = T_a_dist
            P_r_dist_dict[N_slab] = P_r_dist
            P_a_dist_dict[N_slab] = P_a_dist
            x_r_dist_dict[N_slab] = x_r_dist
            Tdp_ao_dist_dict[N_slab] = Tdp_a_dist
            Q_dist_dict[N_slab] = Q_dist
            Q_a_sen_dist_dict[N_slab] = Q_a_sen_dist
            Q_a_lat_dist_dict[N_slab] = Q_a_lat_dist
            SHR_dist_dict[N_slab] = SHR_dist
            DELTAP_air_dist_dict[N_slab] = DELTAP_air_dist
            DELTAP_ref_dist_dict[N_slab] = DELTAP_ref_dist
            htc_ref_dist_dict[N_slab] = htc_ref_dist
            T_ro_tmp = T_ri_tmp
            P_ro_tmp = P_ri_tmp
            x_ro_tmp = x_ri_tmp
            h_ro_tmp = h_ri_tmp
            T_ai_dist[N_slab-1] = np.flip(T_a_dist)
            P_ai_dist[N_slab-1] = np.flip(P_a_dist)
            Tdp_ai_dist[N_slab-1] = np.flip(Tdp_a_dist)
        
        P_ro_guess-=(h_ri_tmp-h_ri)/2e3

    T_ro = T_ro_dict[self.N_slab]; P_ro = P_ro_dict[self.N_slab]; x_ro = x_ro_dict[self.N_slab]
    T_ri = T_ri_tmp; P_ri = P_ri_tmp; x_ri = x_ri_tmp

    T_ao = np.average(T_a_dist_dict[1])
    P_ao = np.average(P_a_dist_dict[1])

    Q_tot = np.sum([np.sum(arr) for arr_list in Q_dist_dict.values() for arr in arr_list])*self.N_tube[0]

    result_tmp = [T_ri,P_ri,x_ri,T_ro,P_ro,x_ro,T_ao,P_ao,Q_tot,T_ro_dict,P_ro_dict,\
                  x_ro_dict,Q_dist_dict,Q_a_sen_dist_dict,Q_a_lat_dist_dict,T_r_dist_dict,\
                  T_a_dist_dict,P_r_dist_dict,P_a_dist_dict,DELTAP_air_dist_dict,\
                  DELTAP_ref_dist_dict,x_r_dist_dict,htc_ref_dist_dict,Tdp_ao_dist_dict,\
                  SHR_dist_dict,T_ai_dist]
    result = result_organizer(result_tmp)
    return  result
import numpy as np
import CoolProp.CoolProp as CP
from scipy import interpolate

class radiator:
    '''
        Unit system of the project: C, kJ, kPa, kW, s, kg ...
        Input Varibles:
            @T_ri               : [C], refrigerant inlet temperature
            @P_ri               : [kPa], refrigerant inlet pressure
            @vol_dor_ref        : [L/min], refrigerant volumetric flow rate
            @T_ai               : [C], air inlet temperature
            @P_ai               : [kPa], air inlet pressure
            @m_dot_air          : [kg/h], air mass flow rate
            @N_segment          : [-], number of segment
    '''
    def __init__(self,T_ri,P_ri,vol_dot_ref,T_ai,P_ai,m_dot_air,N_segment):
        self.fld='INCOMP::MEG[0.5]'
        # -------------- Geometry Info ------------- #
        self.N_slab=1       # number of slabs
        self.N_pass=1       # number of passes
        self.L_tube_pass=700e-3             # [m], tube length in one pass
        self.N_tube_pass=40                 # [-], number of tubes in one pass
        self.t_wall=0.35e-3                 # [m], wall thickness
        self.t_tube=3e-3                    # [m], tube thickness, tube minor
        self.D_tube=35e-3                   # [m], tube depth, tube major
        self.N_port=3                       # [-], number of ports in one tube
        self.Ra_tube=1e-6                   # [m], roughness of tube inner surface
        self.theta_louver=18                # [deg], louver angle
        self.P_louver=1.2e-3                # [m], louver pitch
        self.L_louver=9e-3                  # [m], louver length
        self.N_louverbank=2                 # [-], number of louver sets per fin
        self.h_fin=10e-3                    # [m], fin height
        self.t_fin=0.1e-3                   # [m], fin thickness
        self.P_fin=1.4e-3                   # [m], fin pitch (18 FPI)
        self.N_finptube = 0.0254/18*self.L_tube_pass     # [-], number of fins on each tube
        self.D_fin=35e-3                    # [m], fin depth (same as tube)
        self.Dh_fin = 4*((self.P_fin-self.t_fin+2*(self.h_fin-self.t_fin))*self.D_fin)/((self.h_fin-self.t_fin)*((self.P_fin-self.t_fin)/2)*2)

        self.H_port=self.t_tube-2*self.t_wall       # [m], height of port
        self.W_port=(self.D_tube-(self.N_port+1)*self.t_wall)/self.N_port       # [m], width of port
        self.P_port_tot = 2*(self.H_port+self.W_port)*(self.N_port-2)+2*(self.H_port+2*(self.W_port-self.H_port/2)+np.pi*self.H_port/2)       # [m], total perimeter of ports in one tube
        self.A_port_tot = (self.N_port-2)*(self.H_port*self.W_port)+2*((self.W_port-self.H_port/2)*self.H_port+np.pi*(self.H_port/2)**2/2)      # [m^2], total cross sectional area for refrigeratn in one tube
        self.Dh_port = 4*self.A_port_tot/self.P_port_tot

        self.W_HX = self.L_tube_pass
        self.H_HX = self.t_tube*self.N_tube_pass*self.N_pass+self.h_fin*(self.N_tube_pass-1)*self.N_pass

        self.N_elemptube=N_segment
        self.N_elem = self.N_tube_pass*self.N_pass*N_segment

        self.dx = self.L_tube_pass/self.N_elemptube

        
        # ------------- Area Calculation ----------- #
        self.A_ref_tot = self.P_port_tot*self.L_tube_pass*self.N_tube_pass*self.N_pass
        self.A_ref_elem = self.A_ref_tot/self.N_elem
        self.A_fin_tot = (self.D_fin*self.h_fin)*2*self.N_finptube*self.N_tube_pass*self.N_pass
        self.A_tube_tot = 2*(self.L_tube_pass-self.t_fin*self.N_finptube)*self.D_tube*(self.N_tube_pass-1)*self.N_pass
        self.A_air_tot = self.A_fin_tot+self.A_tube_tot
        self.A_air_elem = self.A_air_tot/self.N_elem
        self.A_front = self.W_HX*self.H_HX      # frontal area of heat exchanger
        self.A_front_cross = self.L_tube_pass*self.t_tube*self.N_tube_pass*self.N_pass + \
                       self.t_fin*self.h_fin*self.N_finptube*self.N_tube_pass*self.N_pass       # frontal cross sectional area of tube+fin
        self.A_air_free = self.A_front-self.A_front_cross           # minimum free flow area

        self.k_w = 155      # [W/m-K]
        self.rho_w = 2710   # [kg/m^3]

        # ------------- Inlet Condition ------------ #
        self.T_ri=T_ri                                                      # [C]
        self.P_ri=P_ri                                                      # [kPa]
        self.vol_dot_ref=vol_dot_ref                                        # [L/min]
        rho_ri=CP.PropsSI('D','T',T_ri+273.15,'P',P_ri*1e3,self.fld)        # [kg/m^3], density of refrigerant
        self.m_dot_ref = self.vol_dot_ref/60*1e-3*rho_ri                    # [kg/s], ref mass flow rate
        self.m_dot_ref_elem = self.m_dot_ref/self.N_tube_pass
        self.G_ref_elem = self.m_dot_ref_elem/self.A_port_tot
        self.vel_ref=self.vol_dot_ref/60*1e-3/(self.N_tube_pass*self.N_pass)/self.A_port_tot        # [m/s], ref flow speed
        self.T_ai=T_ai                                                      # [C]
        self.P_ai=P_ai                                                      # [kPa]
        self.m_dot_air=m_dot_air/(60*60)                                    # [kg/s], air mass flow rate
        self.m_dot_air_elem = self.m_dot_air/self.N_elem
        self.rho_ai = CP.PropsSI('D','T',self.T_ai+273.15,'P',self.P_ai*1e3,'Air')
        self.G_air = self.m_dot_air/self.A_air_free
        self.vel_air = self.m_dot_air/self.rho_ai/self.A_air_free


def airsidehtc_chang_wang(P_air,k_HX,theta_oh,p_fin,p_louver,h_fin,D_fin,L_louver,w_tube,t_tube,t_fin,T_air,Vel_air):
    '''
    {Air side heat transfer coefficient of louver fin microchannel heat exchanger}
    {! Citation: Chang and Wang 1996, A Generalized Heat Transfer Correlation for Louver Fin Geometry}
    
    Input:
        P_air:	Air inlet pressure, P_atm subtract pressure loss upstream
        k_HX:	Conductivity of the tube and fin material,  240 W/m-C for pure aluminum 180 W/m-C for Aluminum alloy [W/m-C]
        theta_oh:	Louver angle of the fins, usually 28 degree for all Al coil tubes [degree]
        p_fin:	Fin pitch [m]
        p_louver:	Louver pitch [m]
        h_fin:	Fin height [m]
        D_fin:	Depth of fin [m]
        L_louver:	Louver length [m]
        w_tube:	Tube depth[m]
        t_tube:	Tube thickness [m]         
        t_fin：	Fin thickness [m]                                                      
        T_air:	Temperature of air used to evaluate the air properties[C] Should we use film temperature here?
        Vel_air:	Air velocity [m/s]
        RH_ai:	Relative humidity at air inlet

    Output: 
        htc:	Heat transfer coefficient [W/m**2-C]
        eta_f:	Fin efficieny [-]
        eta_t:	Total surface efficiency, include unfinned surface, i.e. temperature effectiveness of a finned surface [-] 
    '''
    # {Air property}
    # P_ai = P.P_air - dP

    T_air = T_air + 273.15
    cp_air = CP.PropsSI('C','P',P_air*1e3,'T',T_air,'air')    # [J/kg]    mixture specific heat per unit humid air
    k_air = CP.PropsSI('L','P',P_air*1e3,'T',T_air,'air') # [W/(m-K)]   mixture thermal conductivity 
    mu_air = CP.PropsSI('V','P',P_air*1e3,'T',T_air,'air')   # [kg/(m-s)]    mixture viscosity
    rho_air = CP.PropsSI('D','P',P_air*1e3,'T',T_air,'air') # [kg/m**3]  density of air
    Pr_air = mu_air * cp_air / k_air
    	
    Re_lp = rho_air * Vel_air * p_louver / mu_air
    print(Re_lp)
    p_tube = h_fin + t_tube
    # "! Citation: Chang and Wang 1997, A generalized heat transfer correlation for Iouver fin geometry, Eq.(9)"
    # "Applicable range: 100 < Re_lp < 3000, corrugated fin geometry"
    j = Re_lp**(-0.49) * (theta_oh/90)**0.27 * (p_fin/p_louver)**(-0.14) * (h_fin/p_louver)**(-0.29) * (w_tube/p_louver)**(-0.23) * (L_louver/p_louver)**0.68 * (p_tube/p_louver)**(-0.28) * (t_fin/p_louver)**(-0.05) # "Colburn j-factor"
    # {{"! Citation: Dong et al. 2007, Heat transfer and pressure drop correlations for the multi-louvered fin compact heat exchangers"}
    # {Largely underestimated air side heat transfer coefficient. This paper didn't give applicable Re_lp range. It is developed at 250 < Re_lp < 2600. Face velocity 2 m/s < Vel_face < 18 m/s}
    # j = 0.26712* Re_lp**(-0.1944) * (theta_oh/90)**0.257 * (p_fin/p_louver)**(-0.5177) * (h_fin/p_louver)**(-1.9045) * (w_tube/p_louver)**(-0.2147) * (L_louver/p_louver)**1.7159 * (t_fin/p_louver)**(-0.05) # "Colburn j-factor"}
    St = j * Pr_air**(-2/3) # "Stanton number"
    print(j)
    htc = rho_air * Vel_air * cp_air * St  # "heat transfer coefficient"
    
    # Not sure about the purpose of effectiveness calculated below.    
    # {!m= ( 2 * htc / ( k_HX * t_fin) )**0.5	"M=(h*P)/(k*A)**0.5" {See Incropera and DeWitt page 118 and 122}}
    m = ( 2 * htc / (k_HX * t_fin)  * (1 + t_fin / D_fin) )**0.5   # "Chang and Wang 1996, Eq. (9)" 
    # {!L=H_fin/2}
    L = h_fin / 2 - t_fin    # "Chang and Wang 1996, Eq. (10)"
    eta_f = np.tanh(m*L) / (m*L)	#"fin efficiency" 
    # {!A_f=L*D_fin + t_fin*L*2	"Heat transfer area of one fin [m**2]" }
    A_f = (h_fin/2 ) * D_fin * 2 + t_fin * L * 2	# "Heat transfer area of one fin [m**2]" 
    A_t = A_f + (p_fin - t_fin) * w_tube + p_fin * t_tube	# "Total heat transfer area, including fin and base [m**2]"
    eta_t = 1 - A_f / A_t * (1 - eta_f)	# "Total surface efficiency, i.e. temperature effectiveness "

    return htc, eta_f, eta_t

def airside_dpdz_Chang_Wang_1999(self,T_ai,P_ai,A_front,A_min,L_p,F_p,F_l,F_t,D_h,L_l,F_d,T_p,D_m,theta):
    '''
    !Citation: 
    Yu-Juei Chang, Kuei-Chang Hsu, Yur-Tsai Lin, Chi-Chuan Wang, A generalized friction correlation for louver fin geometry,
    International Journal of Heat and Mass Transfer, Volume 43, Issue 12, 2000, Pages 2237-2243. Yu-Juei Chang, Wen-Jeng
    Chang, Ming-Chia Li, Chi-Chuan Wang, An amendment of the generalized friction correlation for louver fin geometry,
    International Journal of Heat and Mass Transfer, Volume 49, Issues 21–22, 2006, Pages 4250-4253
    
     * Input Variables:
        @rho_air            : air density, [kg/m^3]
        @A_front            : frontal area, [m^2]
        @A_min              : minimum free flow area, [m^2]
        @L_p                : louver pitch, [m]
        @F_p                : fin pitch, [m]
        @F_l                : fin length, [m]
        @F_t                : fin thickness, [m]
        @D_h                : hydralic diameter of fin array, [m]
        @L_l                : louver length, [m]
        @F_d                : fin depth, [m]
        @T_p                : tube thickness, [m]
        @D_m                : major tube diameter, [m]
        @theta              : louver angle, [deg]

     * Output Variables:
        @htc                : heat transfer coefficient, [W/m^2-K]
    '''
    # mu = dynamic viscosity of air
    rho_air = CP.PropsSI('D','T',T_ai+273.15,'P',P_ai*1e3,'Air')
    mu_air = CP.PropsSI('V','T',T_ai+273.15,'P',P_ai*1e3,'Air')
    T_h = T_p-D_m
    Re_Lp = rho_air*self.vel_air*L_p/mu_air
    if Re_Lp<150:
        f1 = 14.39*Re_Lp**(-0.805*F_p/F_l)*(np.log(1.0+(F_p/L_p)))**3.04
        f2 = (np.log((F_t/F_p)**0.48+0.9))**(-1.435)*(D_h/L_p)**(-3.01)*(np.log(0.5*Re_Lp))**(-3.01)
        f3 = (F_p/L_l)**(-0.308)*(F_d/L_l)**(-0.308)*(np.e**(-0.1167*T_p/D_m))*theta**0.35
    elif Re_Lp>=150 and Re_Lp<5e3:
        f1 = 4.97*Re_Lp**(0.6049-1.064/theta**0.2) * (np.log((F_t/F_p)**0.5+0.9))**(-0.527)
        f2 = ((D_h/L_p)*np.log(0.3*Re_Lp))**(-2.966)*(F_p/L_l)**(-0.7931*T_p/T_h)
        f3 = (T_p/D_m)**(-0.0446)*np.log(1.2+(L_p/F_p)**1.4)**(-3.553)*theta**(-0.477)
    f = f1*f2*f3
    return f

def airside_dpdz_Chang_Wang_1996(f,A_c,A,A_front,rho_1,rho_2,rho_l,G_c,K_c,K_e):
    '''
    !Citation:Yu-Juei Chang, Chi-Chuan Wang, Air
              side performance of brazed aluminum
              heat exchangers, Journal of Enhanced
              Heat Transfer, Volume 3, Issue 1,
              1996, Pages 15-28.

         * Input Variables:
            @A_c            : minimum free flow area, [m^2]
            @A              : total surface area, (fin surface + External tube area), [m^2]
            @A_front        : frontal area, (W_HX*H_HX) [m^2]
            @rho_1,rho_2    : air density inlet, outlet, [kg/m^3]
            @rho_l          : louver fin density, [kg/m^3]
            @DELTA_P        : pressure drop, [Pa]
            @G_c            : mass velocity through minimum flow area, [kg/m^2-s]
            @K_c            : abrupt contraction coefficient
            @K_e            : abrupt expansion coefficient

    
    '''
    sigma = A_c/A_front         # contraction ratio of fin array
    rho_m = (rho_1+rho_2)/2
    DELTA_P = (f*A/A_c*rho_1/rho_m + (K_c+1-sigma**2) + 2*(rho_1/rho_2-1) - (1-sigma**2-K_e)*rho_1/rho_2)*G_c**2/(2*rho_l)
    return DELTA_P

def htc_H2OEG(self,Re,T_r,P_r):
    '''
        !Ciatation: from lecture
         * Input Variables:
            @Re         : Reynolds number
            @T_r        : refrigerant temperature, [C]
            @P_r        : refriegrant pressure, [kPa]
        
         * Output Variables:
            @htc        : heat transfer coefficient of H2OEG, [W/m^2-K]
    '''
    if Re<2300:     # laminar for roundtube
        # interpolation data from table
        ratio_ab = np.array([1.0,1.43,2.0,3.0,4.0,8.0])
        Nu_D_q = np.array([3.61,3.73,4.12,4.79,5.33,6.49])
        Nu_D_Ts = np.array([2.98,3.08,3.39,3.96,4.44,5.60])
        Nu_D_ave = (Nu_D_q+Nu_D_Ts)/2
        htc_H2OEG_Laminar=interpolate.UnivariateSpline(ratio_ab,Nu_D_ave,s=0)        # pass all the points forceabley
        
        Nu_D = htc_H2OEG_Laminar(self.W_port/self.H_port)
    else:           # turbulent flow
        # Gnielinski correlation
        f = (0.79*np.log(Re)-1.64)**(-2)
        Pr = CP.PropsSI('PRANDTL','T',T_r+273.15,'P',P_r*1e3,self.fld)
        Nu_D = (f/8)*(Re-1000)*Pr/(1+12.7*(f/8)**(1/2)*(Pr**(2/3)-1))

    k = CP.PropsSI('L','T',T_r+273.15,'P',P_r*1e3,self.fld)     # [W/m-K]
    htc = Nu_D*k/self.Dh_port

    return htc

def dpdz_f_H2OEG(Re,Dh,roughness):
    '''
        Calculate friction factor of H2OEG
        Churchill (1977) correlation
         * Input Variables:
            @Re         : Reynolds number
            @Dh         : hydraulic diameter
            @sigma      : viscosity
        
         * Output variable:
            @dpdz_f     : []
    '''

    A = (-2.456*np.log((7/Re)**0.9+0.27*roughness/Dh))**16; B = (37530/Re)**16
    dpdz_f = 8*((8/Re)**12+1/(A+B)**1.5)**(1/12)

    return dpdz_f

def DP_H2OEG(Re,roughness,rho_ref,G_ref,Vel_i,Vel_o,z_i,z_o,dx,Dh):
    '''
        * Input Variables:
            @Re             : Reynolds number
            @roughness      : roughness of inner tube [m]
            @rho_ref        : density of refrigerant [kg/m^3]
            @G_ref          : mass velocity of refriegrant [kg/m^2-s]
            @Vel_i          : inlet velocity [m/s]
            @Vel_o          : outlet velocity [m/s]
            @z_i            : inlet height [m]
            @z_o            : outlet height [m]
            @dx             : length of each element [m]
            @Dh             : hydraulic diameter, [m]

         * Output Variables
            @DP_total       : [Pa]
    '''
    g=9.81          # [m/s**2]

    # static pressure: friction, (head loss)rho*g*h_L=rho*g*h_L,major + rho*g*h_L_minor
    #   rho*g*h_L_major=f*l/D*rho*V**2/2; rho*g*h_L_minor=K_L*rho*V**2/2
    #   negligible for minor one, only consider the major one
    dpdz_f = dpdz_f_H2OEG(Re,Dh,roughness)
    DP_static = dpdz_f * dx/Dh * G_ref**2/(2*rho_ref)

    # hydrostatic pressure: gravity, rho g h
    DP_hydrostatic = rho_ref*g*(z_i-z_o)
    # dynamic pressure: velocity, 1/2 rho v^2
    DP_dynamic = 1/2*rho_ref*(Vel_i**2-Vel_o**2)
    # total pressure=static + hydrostatic + dynamic
    DP_total = DP_static+DP_hydrostatic+DP_dynamic

    return DP_total

def DP_air(self,T_ai,P_ai,T_ao,P_ao_tmp):
    '''
        Calculate pressure drop through louver fins of air
    '''
    K_c=0.96;K_e=0.6
    rho_ai = CP.PropsSI('D','T',T_ai+273.15,'P',P_ai*1e3,'Air')
    rho_ao = CP.PropsSI('D','T',T_ao+273.15,'P',P_ao_tmp*1e3,'Air')
    f_air = airside_dpdz_Chang_Wang_1999(self,self.T_ai,self.P_ai,self.A_front,self.A_air_free,self.L_louver,self.P_fin,self.h_fin,\
        self.t_fin,self.Dh_fin,self.L_louver,self.D_fin,self.t_tube,self.D_tube,self.theta_louver)
    DELTAP_air = airside_dpdz_Chang_Wang_1996(f_air,self.A_air_free,self.A_air_tot,self.A_front,rho_ai,rho_ao,self.rho_w,self.G_air,K_c,K_e)

    return DELTAP_air

def epsilon_NTU(self,htc_air,htc_ref,T_ri,P_ri,m_dot_ref,T_ai,P_ai,m_dot_air,R_air):
    '''

    '''
    
    cp_ref = CP.PropsSI('C','T',T_ri+273.15,'P',P_ri*1e3,self.fld)/1e3        # [kJ/kg-K]
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
    T_ro = T_ri - Q/C_ref
    T_ao = T_ai + Q/C_air

    return NTU, sigma, Q, T_ro, T_ao

def radiator_tube(self):
    T_r_dist = np.array([self.T_ri])
    T_a_dist = np.array([self.T_ai])
    P_r_dist = np.array([self.P_ri])
    P_a_dist = np.array([self.P_ai])
    Q_dist = np.array([0])
    DELTAP_air_dist = np.array([0])
    DELTAP_H2OEG_dist = np.array([0])
    T_ri = self.T_ri
    P_ri = self.P_ri
    for i in range(self.N_elemptube):
        [T_ro,P_ro,T_ao,P_ao,Q,DELTAP_air,DELTAP_H2OEG] =\
            radiator_tubeelem(self,T_ri,P_ri)
        T_r_dist = np.append(T_r_dist,T_ro)
        T_a_dist = np.append(T_a_dist,T_ao)
        P_r_dist = np.append(P_r_dist,P_ro)
        P_a_dist = np.append(P_a_dist,P_ao)
        Q_dist = np.append(Q_dist,Q)
        DELTAP_air_dist = np.append(DELTAP_air_dist,DELTAP_air)
        DELTAP_H2OEG_dist = np.append(DELTAP_H2OEG_dist,DELTAP_H2OEG)
        T_ri = T_ro
        P_ri = P_ro
    return T_ro,P_ro,T_r_dist,T_a_dist,P_r_dist,P_a_dist,Q_dist,DELTAP_air_dist,DELTAP_H2OEG_dist

def radiator_tubeelem(self,T_ri,P_ri):
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

    rho_ri = CP.PropsSI('D','T',T_ri+273.15,'P',P_ri*1e3,self.fld)
    mu_ri = CP.PropsSI('V','T',T_ri+273.15,'P',P_ri*1e3,self.fld)
    Re_ref = rho_ri*self.vel_ref*self.Dh_port/mu_ri
    htc_ref = htc_H2OEG(self,Re_ref,T_ri,P_ri)/1e3
    [NTU, sigma, Q, T_ro, T_ao] = epsilon_NTU(self,self.htc_air,htc_ref,T_ri,P_ri,\
        self.m_dot_ref_elem,self.T_ai,self.P_ai,self.m_dot_air_elem,self.R_air)

    # Hydraulic of ref and air
    z_i=0; z_o = 0      # neglect hydrostatic pressure change
    DELTAP_H2OEG = DP_H2OEG(Re_ref,self.Ra_tube,rho_ri,self.G_ref_elem,self.vel_ref,self.vel_ref,z_i,z_o,self.dx,self.Dh_port)
    P_ro = P_ri - DELTAP_H2OEG/1e3

    DELTAP_air_tmp = 0
    DELTAP_air = 20
    P_ao_tmp = P_ai - DELTAP_air_tmp
    while abs(DELTAP_air_tmp-DELTAP_air)>0.2:
        DELTAP_air_tmp = DELTAP_air
        DELTAP_air = DP_air(self,self.T_ai,self.P_ai,T_ao,P_ao_tmp)
        P_ao = self.P_ai-DELTAP_air/1e3


    return T_ro, P_ro, T_ao, P_ao, Q, DELTAP_air, DELTAP_H2OEG

def radiator_system(self):
    '''
        Calculate the radiator system
    '''
    # First determine the heat transfer coefficient of air (same for all the elements)
    [htc_air, eta_f, eta_t] = airsidehtc_chang_wang(self.P_ai,self.k_w,self.theta_louver,self.P_fin,self.P_louver,\
        self.h_fin,self.D_fin,self.L_louver,self.D_tube,self.t_tube,self.t_fin,self.T_ai,self.vel_air)
    self.htc_air = htc_air/1e3
    self.R_air = 1 / (self.htc_air * self.A_air_elem * eta_t)

    # Since uniform distribution, only calculate one tube and leads to all the other tubes
    [T_ro,P_ro,T_r_dist,T_a_dist,P_r_dist,P_a_dist,Q_dist,DELTAP_air_dist,DELTAP_H2OEG_dist] = radiator_tube(self)
    Q_tot = self.N_tube_pass*self.N_pass*np.sum(Q_dist)
    T_ao = np.average(T_a_dist)
    P_ao = np.average(P_a_dist)

    return  T_ro,P_ro,T_ao,P_ao,Q_tot




T_ri = 75       # [C]
P_ri = 150      # [kPa]
vol_dot_ref=40  # [L/min]
T_ai = 25       # [C]
P_ai = 99.5     # [kPa]
rh=0.5
m_dot_air = 5000    # [kg/h]
N_segment = 10
self=radiator(T_ri,P_ri,vol_dot_ref,T_ai,P_ai,m_dot_air,N_segment)
[T_ro,P_ro,T_ao,P_ao,Q_tot] = radiator_system(self)
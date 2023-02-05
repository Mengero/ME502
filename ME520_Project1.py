import numpy as np
import CoolProp.CoolProp as CP
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

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
        self.N_finptube = self.L_tube_pass/self.P_fin     # [-], number of fins on each tube
        self.D_fin=35e-3                    # [m], fin depth (same as tube)
        self.Dh_fin = 4*(self.P_fin-self.t_fin)*(self.h_fin-self.t_fin)/(((self.h_fin-self.t_fin)+(self.P_fin-self.t_fin))*2)

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
        self.P_tube = self.h_fin + self.t_tube
        
        # ------------- Area Calculation ----------- #
        self.A_ref_tot = self.P_port_tot*self.L_tube_pass*self.N_tube_pass*self.N_pass
        self.A_ref_elem = self.A_ref_tot/self.N_elem
        self.A_fin_tot = (self.D_fin*self.h_fin)*2*self.N_finptube*self.N_tube_pass*self.N_pass
        self.A_tube_tot = 2*(self.L_tube_pass-self.t_fin*self.N_finptube)*self.D_tube*(self.N_tube_pass-1)*self.N_pass
        self.A_air_tot = self.A_fin_tot+self.A_tube_tot
        self.A_air_elem = self.A_air_tot/self.N_elem
        self.A_front = self.W_HX*self.H_HX      # frontal area of heat exchanger
        self.A_front_cross = self.L_tube_pass*self.t_tube*self.N_tube_pass*self.N_pass + \
                       (self.t_fin*self.h_fin+self.t_fin*(self.P_fin-self.t_fin))*self.N_finptube*self.N_tube_pass*self.N_pass # frontal cross sectional area of tube+fin
        self.A_air_free = self.A_front-self.A_front_cross           # minimum free flow area
        self.sigma = self.A_air_free/self.A_front
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
        # self.vel_air = 4.758


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
    Pr_air = CP.PropsSI('PRANDTL','P',P_air*1e3,'T',T_air,'air')
    	
    Re_lp = rho_air * Vel_air * p_louver / mu_air
    p_tube = h_fin + t_tube
    # "! Citation: Chang and Wang 1997, A generalized heat transfer correlation for Iouver fin geometry, Eq.(9)"
    # "Applicable range: 100 < Re_lp < 3000, corrugated fin geometry"
    j = Re_lp**(-0.49) * (theta_oh/90)**0.27 * (p_fin/p_louver)**(-0.14) * (h_fin/p_louver)**(-0.29) * (w_tube/p_louver)**(-0.23) * (L_louver/p_louver)**0.68 * (p_tube/p_louver)**(-0.28) * (t_fin/p_louver)**(-0.05) # "Colburn j-factor"
    # {{"! Citation: Dong et al. 2007, Heat transfer and pressure drop correlations for the multi-louvered fin compact heat exchangers"}
    # {Largely underestimated air side heat transfer coefficient. This paper didn't give applicable Re_lp range. It is developed at 250 < Re_lp < 2600. Face velocity 2 m/s < Vel_face < 18 m/s}
    # j = 0.26712* Re_lp**(-0.1944) * (theta_oh/90)**0.257 * (p_fin/p_louver)**(-0.5177) * (h_fin/p_louver)**(-1.9045) * (w_tube/p_louver)**(-0.2147) * (L_louver/p_louver)**1.7159 * (t_fin/p_louver)**(-0.05) # "Colburn j-factor"}
    St = j * Pr_air**(-2/3) # "Stanton number"
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

def airside_dpdz_Chang_Wang_1999(self,T_ai,P_ai,L_p,F_p,F_l,F_t,D_h,L_l,F_d,T_p,D_m,theta):
    '''
    !Citation: 
    Yu-Juei Chang, Kuei-Chang Hsu, Yur-Tsai Lin, Chi-Chuan Wang, A generalized friction correlation for louver fin geometry,
    International Journal of Heat and Mass Transfer, Volume 43, Issue 12, 2000, Pages 2237-2243. Yu-Juei Chang, Wen-Jeng
    Chang, Ming-Chia Li, Chi-Chuan Wang, An amendment of the generalized friction correlation for louver fin geometry,
    International Journal of Heat and Mass Transfer, Volume 49, Issues 21–22, 2006, Pages 4250-4253
    
     * Input Variables:
        @T_ai, P_ai
        @L_p                : louver pitch, [m]
        @F_p                : fin pitch, [m]
        @F_l                : fin length, [m]
        @F_t                : fin thickness, [m]
        @D_h                : hydralic diameter of fin array, [m]
        @L_l                : louver length, [m]
        @F_d                : fin depth, [m]
        @T_p                : tube pitch, [m]
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
        f2 = (((D_h/L_p)*np.log(0.3*Re_Lp))**(-2.966)) * (F_p/L_l)**(-0.7931*T_p/T_h)
        f3 = ((T_p/D_m)**(-0.0446)) * ((np.log(1.2+(L_p/F_p)**1.4))**(-3.553)) * theta**(-0.477)
    f = f1*f2*f3
    return f

def airside_dpdz_Chang_Wang_1996(f,A_c,A,A_front,rho_1,rho_2,G_c,K_c,K_e):
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
    DELTA_P = (f*A/A_c*rho_1/rho_m + (K_c+1-sigma**2) + 2*(rho_1/rho_2-1) - (1-sigma**2-K_e)*rho_1/rho_2)*G_c**2/(2*rho_1)
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
    sigma_K_c=np.array([0.0998923586380184,0.1533285509902264,0.20114352910644695,0.2559862929607125,0.30027807995879835,\
        0.35159596210418154,0.4008039869964784,0.4528217737113951,0.5006147741478959,0.5540272982296366,0.6018135363031468,\
        0.6481915122838517,0.7030004643231641,0.7521797491727505,0.8027672961151423,0.856154461335668,0.902517221999644,\
        0.9488799826636198,0.9980457427872251])
    K_c_tmp=np.array([0.39802350427350386,0.3896768162393158,0.3844601362179483,0.37611344818376025,0.3667234241452988,\
        0.3531600560897432,0.3395966880341874,0.323946647970085,0.3051665998931621,0.2822132077991448,0.25925981570512757,\
        0.23526308760683734,0.20604967948717912,0.17474959935897405,0.1444928552350424,0.10588942307692273,0.0725026709401706,\
        0.03911591880341825,-0.0005308493589750718])
    sigma_K_e=np.array([0.10196840407615151,0.14404551719479566,0.20156279561186796,0.24785962323668492,0.2997776650974894,\
        0.3531090408233421,0.4008344176299364,0.4499731283015789,0.5012250774078031,0.5496790988266411,0.6016529301821186,\
        0.653631833309839,0.7028128087501732,0.7527021417137791,0.8025982370403758,0.8532060710717395,0.9024124053732887,\
        0.9488157402152091,0.9987473379475063])
    K_e_tmp=np.array([0.8111845619658116,0.7371077056623927,0.6432074652777773,0.5691306089743584,0.49192374465811906,\
        0.41889022435897405,0.35837673611111076,0.3020365918803414,0.24778311965811906,0.2029196714743584,0.16014289529914527,\
        0.12049612713675151,0.09023938301282008,0.06311264690170915,0.04015925480769189,0.022422542735042184,0.007815838675213183,\
        -0.0005308493589748497,-0.0015741853632484926])

    K_c_interpolate=interpolate.UnivariateSpline(sigma_K_c,K_c_tmp,s=0)        # pass all the points forceabley
    K_e_interpolate=interpolate.UnivariateSpline(sigma_K_e,K_e_tmp,s=0)
    K_c=K_c_interpolate(self.sigma)
    K_e=K_e_interpolate(self.sigma)
    rho_ai = CP.PropsSI('D','T',T_ai+273.15,'P',P_ai*1e3,'Air')
    rho_ao = CP.PropsSI('D','T',T_ao+273.15,'P',P_ao_tmp*1e3,'Air')
    f_air = airside_dpdz_Chang_Wang_1999(self,self.T_ai,self.P_ai,self.P_louver,self.P_fin,self.h_fin,\
        self.t_fin,self.Dh_fin,self.L_louver,self.D_fin,self.P_tube,self.t_tube,self.theta_louver)
    DELTAP_air = airside_dpdz_Chang_Wang_1996(f_air,self.A_air_free,self.A_air_tot,self.A_front,rho_ai,rho_ao,self.G_air,K_c,K_e)

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
    P_r_dist = np.array([self.P_ri])
    T_a_dist = np.array([])
    P_a_dist = np.array([])
    Q_dist = np.array([])
    DELTAP_air_dist = np.array([])
    DELTAP_H2OEG_dist = np.array([])
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
    P_ao_tmp = self.P_ai - DELTAP_air_tmp
    while abs(DELTAP_air_tmp-DELTAP_air)>1e-8:
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

    return  T_ro,P_ro,T_ao,P_ao,Q_dist,Q_tot,T_r_dist,T_a_dist,P_r_dist,P_a_dist,DELTAP_air_dist,DELTAP_H2OEG_dist


def N_segment_effect(self):
    T_ri = 75       # [C]
    P_ri = 150      # [kPa]
    vol_dot_ref=40  # [L/min]
    T_ai = 25       # [C]
    P_ai = 99.5     # [kPa]
    rh=0.5
    m_dot_air = 5000    # [kg/h]

    # Define data storage
    N_segment_1_dict = {}
    N_segment_2_dict = {}
    T_ro_dict = {}
    P_ro_dict = {}
    T_ao_dict = {}
    P_ao_dict = {}
    Q_dist_dict = {}
    Q_dict = {}
    T_r_dist_dict = {}
    T_a_dist_dict = {}
    P_r_dist_dict = {}
    P_a_dist_dict = {}
    DELTAP_air_dist_dict = {}
    DELTAP_H2OEG_dist_dict = {}

    for N_segment in range(1,51):
        N_segment_1 = np.array(range(N_segment+1))
        N_segment_2 = np.array(range(1,N_segment+1))
        self=radiator(T_ri,P_ri,vol_dot_ref,T_ai,P_ai,m_dot_air,N_segment)
        [T_ro,P_ro,T_ao,P_ao,Q_dist,Q_tot,T_r_dist,T_a_dist,P_r_dist,P_a_dist,DELTAP_air_dist,DELTAP_H2OEG_dist] = radiator_system(self)
        N_segment_1_dict[N_segment] = N_segment_1
        N_segment_2_dict[N_segment] = N_segment_2
        T_ro_dict[N_segment] = T_ro; P_ro_dict[N_segment] = P_ro
        T_ao_dict[N_segment] = T_ao; P_ao_dict[N_segment] = P_ao
        Q_dist_dict[N_segment] = Q_dist
        Q_dict[N_segment] = Q_tot
        T_r_dist_dict[N_segment] = T_r_dist; T_a_dist_dict[N_segment] = T_a_dist
        P_r_dist_dict[N_segment] = P_r_dist; P_a_dist_dict[N_segment] = P_a_dist
        DELTAP_air_dist_dict[N_segment] = DELTAP_air_dist
        DELTAP_H2OEG_dist_dict[N_segment] = DELTAP_H2OEG_dist

    return N_segment_1_dict,N_segment_2_dict,T_ro_dict,P_ro_dict,T_ao_dict,P_ao_dict,Q_dist_dict,Q_dict,T_r_dist_dict,T_a_dist_dict,P_r_dist_dict,P_a_dist_dict,DELTAP_air_dist_dict,DELTAP_H2OEG_dist_dict

def plot_result_segment(self):
    [N_segment_1_dict,N_segment_2_dict,T_ro_dict,P_ro_dict,T_ao_dict,P_ao_dict,Q_dist_dict,Q_dict,T_r_dist_dict,\
        T_a_dist_dict,P_r_dist_dict,P_a_dist_dict,DELTAP_air_dist_dict,DELTAP_H2OEG_dist_dict] = N_segment_effect(self)
    print(T_ro_dict[40],P_ro_dict[40],T_ao_dict[40],P_ao_dict[40],Q_dict[40])
    # ------------ plot overall performance of the heat exchanger at 50 segments ------------- #
    N_elemptube = 40
    x_1 = N_segment_1_dict[N_elemptube]; x_2 = N_segment_2_dict[N_elemptube]
    overall_performance_1 = [T_r_dist_dict[N_elemptube],P_r_dist_dict[N_elemptube]]
    overall_performance_2 = [T_a_dist_dict[N_elemptube],P_a_dist_dict[N_elemptube],\
                             DELTAP_air_dist_dict[N_elemptube],\
                             Q_dist_dict[N_elemptube],DELTAP_H2OEG_dist_dict[N_elemptube]]
    
    fig1, (ax1,ax2) = plt.subplots(2,1)
    ax1.plot(x_1,overall_performance_1[0],'o-')
    ax1.grid()
    ax1.set_title('Radiator H2O-EG Temp Distribution, N=40')
    ax1.set_ylabel('T [C]')
    ax2.plot(x_1,overall_performance_1[1],'o-')
    ax2.grid()
    ax2.set_title('Radiator H2O-G Pressure Distribution, N=40')
    ax2.set_xlabel('Location [-]')
    ax2.set_ylabel('P [kPa]')
    fig1.tight_layout()
    fig1.set_figheight(14)
    fig1.set_figwidth(14)
    fig1.savefig('Radiator_Overall_Performance_1.pdf')

    fig2, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,1)
    ax1.plot(x_2,overall_performance_2[0],'o-')
    ax1.grid()
    ax1.set_title('Radiator Air Temperature Outlet Distribution, N=40')
    ax1.set_ylabel('T [C]')
    ax2.plot(x_2,overall_performance_2[1],'o-')
    ax2.grid()
    ax2.set_title('Radiator Air Pressure Outlet Distribution, N=40')
    ax2.set_ylabel('P [kPa]')
    ax3.plot(x_2,overall_performance_2[2],'o-')
    ax3.grid()
    ax3.set_title('Radiator Air Pressure Drop Distribution, N=40')
    ax3.set_ylabel('P [Pa]')
    ax4.plot(x_2,overall_performance_2[3],'o-')
    ax4.grid()
    ax4.set_title('Radiator Heat Transfer Rate Distribution, N=40')
    ax4.set_ylabel('Q [kW]')
    ax5.plot(x_2,overall_performance_2[4],'o-')
    ax5.grid()
    ax5.set_title('Radiator H2O-EG Pressure Drop Distribution, N=40')
    ax5.set_ylabel('P [Pa]')
    ax5.set_xlabel('Location (flow from left to right) [-]')
    fig2.set_figheight(18)
    fig2.set_figwidth(14)
    fig2.savefig('Radiator_Overall_Performance_2.pdf')

    # ----------- plot segment impact on the model ----------- #
    fig3, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,1)
    T_ro_list=T_ro_dict.items()
    x,y = zip(*T_ro_list)
    ax1.plot(x,y,'o-')
    ax1.grid()
    ax1.set_title('Radiator H2O-EG Temperature Outlet')
    ax1.set_ylabel('T [C]')

    P_ro_list=P_ro_dict.items()
    x,y=zip(*P_ro_list)
    ax2.plot(x,y,'o-')
    ax2.grid()
    ax2.set_title('Radiator H2O-EG Pressure Outlet')
    ax2.set_ylabel('P [kPa]')

    T_ao_list=T_ao_dict.items()
    x,y=zip(*T_ao_list)
    ax3.plot(x,y,'o-')
    ax3.grid()
    ax3.set_title('Radiator Air Temperature Outlet')
    ax3.set_ylabel('T [C]')

    P_ao_list=P_ao_dict.items()
    x,y=zip(*P_ao_list)
    ax4.plot(x,y,'o-')
    ax4.grid()
    ax4.set_title('Radiator Air Pressure Outlet')
    ax4.set_ylabel('P [kPa]')

    Q_list=Q_dict.items()
    x,y=zip(*Q_list)
    ax5.plot(x,y,'o-')
    ax5.grid()
    ax5.set_title('Radiator Heat Transfer Rate')
    ax5.set_ylabel('Q [kW]')
    ax5.set_xlabel('N_segment [-]')
    fig3.set_figheight(18)
    fig3.set_figwidth(14)
    fig3.savefig('Radiator_Segment_Impact.pdf')

    return

def plot_result_m_dot_air(self):
    T_ri = 75       # [C]
    P_ri = 150      # [kPa]
    vol_dot_ref=40  # [L/min]
    T_ai = 25       # [C]
    P_ai = 99.5     # [kPa]
    N_segment = 40
    T_ro_dict={}
    P_ro_dict={}
    T_ao_dict={}
    P_ao_dict={}
    Q_tot_dict={}
    for m_dot_air in range(4000,6000,100):
        self=radiator(T_ri,P_ri,vol_dot_ref,T_ai,P_ai,m_dot_air,N_segment)
        result_tmp = radiator_system(self)
        # T_ro,P_ro,T_ao,P_ao,Q_dist,Q_tot,T_r_dist,T_a_dist,P_r_dist,P_a_dist,DELTAP_air_dist,DELTAP_H2OEG_dist
        T_ro_dict[m_dot_air]=result_tmp[0]
        P_ro_dict[m_dot_air]=result_tmp[1]
        T_ao_dict[m_dot_air]=result_tmp[2]
        P_ao_dict[m_dot_air]=result_tmp[3]
        Q_tot_dict[m_dot_air]=result_tmp[5]
    fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,1)
    T_ro_list=T_ro_dict.items()
    x,y=zip(*T_ro_list)
    ax1.plot(x,y,'-')
    ax1.set_ylabel('T [C]')
    ax1.set_title('Outlet coolant temperature vs. air mass flow rate')
    
    P_ro_list=P_ro_dict.items()
    x,y=zip(*P_ro_list)
    ax2.plot(x,y,'-')
    ax2.set_ylabel('P [kPa]')
    ax2.set_title('Outlet coolant pressure vs. air mass flow rate')

    T_ao_list=T_ao_dict.items()
    x,y=zip(*T_ao_list)
    ax3.plot(x,y,'-')
    ax3.set_ylabel('T [C]')
    ax3.set_title('Outlet air emperature vs. air mass flow rate')

    P_ao_list=P_ao_dict.items()
    x,y=zip(*P_ao_list)
    ax4.plot(x,y,'-')
    ax4.set_ylabel('P [kPa]')
    ax4.set_title('Outlet air pressure vs. air mass flow rate')

    Q_tot_list=Q_tot_dict.items()
    x,y=zip(*Q_tot_list)
    ax5.plot(x,y,'-')
    ax5.set_ylabel('Q [kW]')
    ax5.set_title('Heat transfer rate vs. air mass flow rate')
    ax5.set_xlabel('air mass flow rate [kg/h]')
    fig.set_figheight(18)
    fig.set_figwidth(14)
    fig.savefig('Radiator_airmassflowrate_Impact.pdf')
    return

T_ri = 75       # [C]
P_ri = 150      # [kPa]
vol_dot_ref=40  # [L/min]
T_ai = 25       # [C]
P_ai = 99.5     # [kPa]
rh=0.5
m_dot_air = 5000    # [kg/h]
N_segment = 40
self=radiator(T_ri,P_ri,vol_dot_ref,T_ai,P_ai,m_dot_air,N_segment)
# [T_ro,P_ro,T_ao,P_ao,Q_tot,T_r_dist,T_a_dist,P_r_dist,P_a_dist,DELTAP_air_dist,DELTAP_H2OEG_dist] = radiator_system(self)

plot_result_segment(self)

plot_result_m_dot_air(self)
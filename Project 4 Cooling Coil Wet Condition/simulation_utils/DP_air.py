import numpy as np
from scipy import interpolate
import CoolProp.CoolProp as CP

def airside_dpdz_Chang_Wang_1999(vel_air,T_ai,P_ai,L_p,F_p,F_l,F_t,D_h,L_l,F_d,T_p,D_m,theta):
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
    Re_Lp = rho_air*vel_air*L_p/mu_air
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

def DP_air(sigma,T_ai,P_ai,T_ao,P_ao_tmp,vel_air,P_louver,P_fin,h_fin,t_fin,Dh_fin,L_louver,D_fin,\
           P_tube,t_tube,theta_louver,A_air_free,A_air_tot,A_front,G_air):
    '''
        Calculate pressure drop through louver fins of air
        
         * Input Variables:
            @sigma                      : contraction ratio of air flow area, [-]
            @T_ai,P_ai,T_ao,P_ao_tmp    : air inlet & outlet properties, [C] [kPa]
            @vel_air                    : air velocity (contracted), [m/s]
            @P_louver                   : louver pitch, [m]
            @P_fin                      : fin pitch, [m]
            @h_fin                      : fin height, [m]
            @t_fin                      : fin thickness, [m]
            @Dh_fin                     : hydraulic diameter, [m]
            @L_louver                   : louver length, [m]
            @D_fin                      : fin depth, [m]
            @P_tube                     : tube pitch, [m]
            @t_tube                     : tube thickness, [m]
            @theta_louver               : louver angle, [rad]
            @A_air_free                 : air side free flow area, [m^2]

         * Output Variables:
            @DELTAP_air                 : air side pressure drop, [Pa]
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
    K_c=K_c_interpolate(sigma)
    K_e=K_e_interpolate(sigma)
    rho_ai = CP.PropsSI('D','T',T_ai+273.15,'P',P_ai*1e3,'Air')
    rho_ao = CP.PropsSI('D','T',T_ao+273.15,'P',P_ao_tmp*1e3,'Air')
    f_air = airside_dpdz_Chang_Wang_1999(vel_air,T_ai,P_ai,P_louver,P_fin,h_fin,\
        t_fin,Dh_fin,L_louver,D_fin,P_tube,t_tube,theta_louver)
    DELTAP_air = airside_dpdz_Chang_Wang_1996(f_air,A_air_free,A_air_tot,A_front,rho_ai,rho_ao,G_air,K_c,K_e)

    return DELTAP_air
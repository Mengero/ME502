import numpy as np
import CoolProp.CoolProp as CP

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
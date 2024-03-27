import numpy as np
import CoolProp.CoolProp as CP

def dpdz_f_2ph_Del_Col_2013(fld,roughness,Dh,G_ref,P_r,x_r):
    '''
        !Citation: Del Col, Davide, et al. "Experiments and updated model for two phase frictional pressure drop inside minichannels." International journal of heat and mass transfer 67 (2013): 326-337.
         * Input Variables
            @fld            : fluid name
            @roughness      : roughness of inner tube wall, [m]
            @Dh             : tube inner hydraulic diameter, [m]
            @G_ref          : mass flux, [kg/m^2-s]
            @P_r            : refrigerant pressure, [kPa]
            @x_r            : refrigerant quality, [-]

         * Output Variables
            @dpdz_f         : [Pa/m]
    
    '''
    mu_l = CP.PropsSI('V','P',P_r*1e3,'Q',0,fld)          # N/m^2-s
    P_crit = CP.PropsSI('PCRIT','P',P_r*1e3,'Q',x_r,fld)/1e3        # [kPa]
    rho_liq = CP.PropsSI('D','P',P_r*1e3,'Q',0,fld)
    rho_gas = CP.PropsSI('D','P',P_r*1e3,'Q',1,fld)
    mu_liq = CP.PropsSI('V','P',P_r*1e3,'Q',0,fld)
    mu_gas = CP.PropsSI('V','P',P_r*1e3,'Q',1,fld)

    RR = 2*roughness/Dh
    Re_LO = G_ref * Dh / mu_l
    A=0.046*(3500**(-0.2))
    Re_LO_plus=((A+0.7*RR)/0.046)**(-5)

    if Re_LO<=Re_LO_plus:
        X = 0
    elif Re_LO>=3500:
        X=1
    elif Re_LO<3500 and Re_LO>Re_LO_plus:
        X=1+(A-0.046*Re_LO**(-0.2))/(0.7*RR)

    Z = (1-x_r)**2+x_r**2*rho_liq/rho_gas*(mu_gas/mu_liq)**(0.2)
    W=1.398*(P_r/P_crit)
    F = x_r**0.9525*(1-x_r)**0.414
    H = (rho_liq/rho_gas)**1.132 * (mu_gas/mu_liq)**0.44 * (1-mu_gas/mu_liq)**3.542
    rho_mean = rho_gas * (1 + G_ref*x_r/(Dh**2*np.pi/4*rho_gas*(G_ref*x_r/rho_gas)))
    E=0.015+0.44*np.log((rho_mean/rho_liq)*(mu_liq))
    phi_LO = np.sqrt(Z+3.595*F*H*(1-E)**W)
    f_LO = 0.046*(Re_LO)**(-0.2)+0.7*RR*X
    dpdz_f = phi_LO**2 * (2*f_LO*G_ref**2/(Dh*rho_liq))
    return  f_LO
import numpy as np
import CoolProp.CoolProp as CP

def friction_factor(Re_k):
    '''
        Calculate friction factor based on Reynolds number
    '''
    
    if Re_k<2e3:
        f_k = 16*Re_k**(-1)
    elif Re_k>=2e3 and Re_k<2e4:
        f_k = 0.079*Re_k**(-0.25)
    elif Re_k>2e4:
        f_k = 0.046*Re_k**(-0.2)
    
    return f_k

def htc_2ph_r1234yf_Kim_2014(fld,G_ref,Dh_port,asp_ratio,P_r,x_r):
    '''
        !Citation: Kim, Sung-Min, and Issam Mudawar. "Review of databases and predictive methods for heat transfer in condensing and boiling mini/micro-channel flows." International Journal of Heat and Mass Transfer 77 (2014): 627-652.

         * Input Variables
            @self
            @P_r            : refrigerant pressure, [kPa]
            @x_r            : refrigearnt quality, [0-1]

         * Output Variables
            @htc_ref        : refrigerant heat transfer coefficient, [W/m^2-K]

    '''
    mu_f = CP.PropsSI('V','P',P_r*1e3,'Q',0,fld)       # [N/m^2-s]
    mu_g = CP.PropsSI('V','P',P_r*1e3,'Q',1,fld)
    rho_f = CP.PropsSI('D','P',P_r*1e3,'Q',0,fld)      # [kg/m^3]
    rho_g = CP.PropsSI('D','P',P_r*1e3,'Q',1,fld)
    sigma = CP.PropsSI('I','P',P_r*1e3,'Q',x_r,fld)    # [N/m]
    Pr_f = CP.PropsSI('PRANDTL','P',P_r*1e3,'Q',0,fld)
    k_f = CP.PropsSI('L','P',P_r*1e3,'Q',0,fld)        # [W/m-K]
    X_tt = (mu_f/mu_g)**0.1 * ((1-x_r)/x_r)**0.9 * (rho_g/rho_f)**0.5
    Re_f = G_ref*(1-x_r)*Dh_port/mu_f
    Re_g = G_ref*x_r*Dh_port/mu_g
    Re_fo = G_ref*Dh_port/mu_f
    Su_go = rho_g*sigma*Dh_port/mu_g**2

    if Re_f >=2000 and Re_g>=2000:
        C=0.39*Re_fo**0.03*Su_go**0.1*(rho_f/rho_g)**0.35
    elif Re_f>=2000 and Re_g<2000:
        C=8.7e-4*Re_fo**0.17*Su_go**0.5*(rho_f/rho_g)**0.36
    elif Re_f<2000 and Re_g>=2000:
        C=0.0015*Re_fo**0.59*Su_go**0.19*(rho_f/rho_g)**0.36
    elif Re_f<2000 and Re_g<2000:
        C=3.5e-5*Re_fo**0.44*Su_go**0.5*(rho_f/rho_g)**(0.48)

    f_f = friction_factor(Re_f); f_g = friction_factor(Re_g)
    v_f = mu_f/rho_f; v_g = mu_g/rho_g
    dpdz_f = -2*f_f*v_f*G_ref**2*(1-x_r)**2/Dh_port
    dpdz_g = -2*f_g*v_g*G_ref**2*x_r**2/Dh_port
    X = np.sqrt(dpdz_f/dpdz_g)
    phi_g = np.sqrt(1+C*X+X**2)
    
    if Re_f<=1250:
        We=2.45*Re_g**0.64/(Su_go**0.3*(1+1.09*X_tt**0.039))**0.4
    else:
        We=0.85*(Re_g**0.79*X_tt**0.157)/(Su_go**0.3*(1+1.09*X_tt**0.039))**0.4*((mu_g/mu_f)**2*(v_g/v_f))**0.084

    if We>7*X_tt**0.2:
        htp_cir=(0.048*Re_f**0.69*Pr_f**0.34*phi_g/X_tt)*k_f/Dh_port
    else:
        htp_cir=((0.048*Re_f**0.69*Pr_f**0.34*phi_g/X_tt)**2+(3.2e-7*Re_f**(-0.38)*Su_go**1.39))

    Nu_3 = 8.235*(1-1.833*asp_ratio+3.767*asp_ratio**2-5.814*asp_ratio**3+5.361*asp_ratio**4-2*asp_ratio**5)
    Nu_4 = 8.235*(1-2.042*asp_ratio+3.085*asp_ratio**2-2.477*asp_ratio**3+1.058*asp_ratio**4-0.186*asp_ratio**5)
    htc_ref = htp_cir*(Nu_3/Nu_4)

    return htc_ref
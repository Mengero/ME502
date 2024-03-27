import numpy as np  
import CoolProp.CoolProp as CP

def f_phase(Re):
    '''
        calculate liquid/gas only friction factor
    '''
    if Re<=2000:
        f = 64/Re
    elif Re>=3000:
        f = 0.25*(np.log(150.39/Re**0.98865 - 152.66/Re))**(-2)
    elif Re>2000 and Re<3000:
        f = (1.1525*Re+895)*1e-5

    return f

def dpdz_f_2ph_Lopez_2014(fld,G_ref,T_sat,x_r,Dh):
    '''
        !Citation: Lopez-Belchi, Alejandro, et al. "Experimental condensing two-phase frictional pressure drop inside mini-channels. Comparisons and new model development." International Journal of Heat and Mass Transfer 75 (2014): 581-591.
         * Input Variables
            @fld                : fluid name
            @G_ref              : mass flux, [kg/m^2-s]
            @P_r                : refrigerant pressure, [kPa]
            @x_r                : refrigerant quality, [-]
            @Dh                 : hydraulic diameter, [m]

         * Output Variables
            @dpdz_f             : [Pa/m]
    '''
    rho_liq = CP.PropsSI('D','T',T_sat+273.15,'Q',0,fld)
    rho_gas = CP.PropsSI('D','T',T_sat+273.15,'Q',1,fld)
    mu_liq = CP.PropsSI('V','T',T_sat+273.15,'Q',0,fld)       # [N/m^2-s]
    mu_gas = CP.PropsSI('V','T',T_sat+273.15,'Q',1,fld)
    P_crit = CP.PropsSI('PCRIT','T',T_sat+273.15,'Q',x_r,fld)/1e3
    P_r = CP.PropsSI('P','T',T_sat+273.15,'Q',x_r,fld)/1e3

    Re_liq = G_ref*(1-x_r)*Dh/mu_liq
    Re_gas = G_ref*x_r*Dh/mu_gas
    f_liq = f_phase(Re_liq); f_gas = f_phase(Re_gas)

    dpdz_liq = G_ref**2 * (1-x_r)**2 / (2 * Dh * rho_liq) * f_liq
    dpdz_gas = G_ref**2 * x_r**2 / (2 * Dh * rho_gas) * f_gas 
    X = np.sqrt(dpdz_liq / dpdz_gas)
    C=4.6468e-6*(P_r/P_crit)**5.5866*Re_liq**0.4387*(rho_liq/rho_gas)**5.7189*X**(-0.4243)
    phi_liq_2 = 1+C/X+1/X**2
    dpdz_f = phi_liq_2 * dpdz_liq

    return dpdz_f

# T_sat = 35
# fld='R1234yf'
# G_ref=475       # [kg/m^2-s]
# x_r = 0.4
# Dh = 1.16e-3
# P_r = CP.PropsSI('P','T',T_sat+273.15,'Q',x_r,fld)/1e3

# dpdz = dpdz_f_2ph_Lopez_2014(fld,G_ref,T_sat,x_r,Dh)
# print(dpdz)
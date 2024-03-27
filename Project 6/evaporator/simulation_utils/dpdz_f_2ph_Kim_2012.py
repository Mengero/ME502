import numpy as np
import CoolProp.CoolProp as CP

def f_phase(Re):
    '''
        calculate liquid/gas only friction factor
    '''
    if Re<2000:
        f = 16/Re
    elif Re>=20000:
        f = 0.046*Re**(-0.2)
    elif Re>=2000 and Re<20000:
        f = 0.079*Re**(-0.25)

    return f

def dpdz_f_2ph_Kim_2012(Ref,G,x,T_sat,Dh,q_H,P_H,P_F):
    '''
        !Citation: Kim, Sung-Min, and Issam Mudawar. "Universal approach to predicting two-phase frictional pressure drop for mini/micro-channel saturated flow boiling." International Journal of Heat and Mass Transfer 58.1-2 (2013): 718-734.
         * Input Variable
            @Ref            : refrigerant name
            @G              : refrigerant mass flux, [kg/m^2-s]
            @x              : refrigerant quality, [-]
            @T_sat          : refrigerant satuartion temperature, [C]
            @Dh             : hydraulic diameter, [m]
            @q_H            : effective heat flux, [W/m^2]
            @P_H            : heated perimeter of channel
            @P_F            : wetted perimeter of channel

         * Output Variable
            @dpdz_f         : Fanning friction factor, [-]
    '''
    mu_f = CP.PropsSI('V','T',T_sat+273.15,'Q',0,Ref)
    mu_g = CP.PropsSI('V','T',T_sat+273.15,'Q',1,Ref)
    rho_f = CP.PropsSI('D','T',T_sat+273.15,'Q',0,Ref)
    rho_g = CP.PropsSI('D','T',T_sat+273.15,'Q',1,Ref)
    sigma = CP.PropsSI('I','T',T_sat+273.15,'Q',x,Ref)
    h_f = CP.PropsSI('H','T',T_sat+273.15,'Q',0,Ref)
    h_g = CP.PropsSI('H','T',T_sat+273.15,'Q',1,Ref)
    h_fg = h_g - h_f
    Re_f = G*(1-x)*Dh/mu_f; Re_g = G*x*Dh/mu_g
    Re_fo = G*Dh/mu_f
    Su_go = rho_g*sigma*Dh/mu_g**2
    We_fo = G**2*Dh/(rho_f*sigma)
    Bo=q_H/(G*h_fg)
    if Re_f>=2000 and Re_g>=2000:               # tt
        C_nonboiling = 0.39*Re_fo**0.03*Su_go**0.1*(rho_f/rho_g)**0.35
    elif Re_f>=2000 and Re_g<2000:              # tv
        C_nonboiling = 8.7e-4*Re_fo**0.59*Su_go**0.5*(rho_f/rho_g)**0.14
    elif Re_f<2000 and Re_g>=2000:
        C_nonboiling = 1.5e-3*Re_fo**0.59*Su_go**0.19*(rho_f/rho_g)**0.36
    elif Re_f<2000 and Re_g<2000:
        C_nonboiling = 3.5e-5*Re_fo**0.44*Su_go**0.5*(rho_f/rho_g)**0.48

    if Re_f>=2000:
        C = C_nonboiling*(1+60*We_fo**0.32*(Bo*P_H/P_F)**0.78)
    else:
        C = C_nonboiling*(1+530*We_fo**0.52*(Bo*P_H/P_F)**1.09)

    f_f = f_phase(Re_f); f_g = f_phase(Re_g)
    dpdz_f = -(2*f_f*G**2*(1-x)**2)/(Dh*rho_f)
    dpdz_g = -(2*f_g*G**2*x**2)/(Dh*rho_g)
    X = np.sqrt(dpdz_f/dpdz_g)
    phi_f_2 = 1+C/X+1/X**2
    dpdz_F = dpdz_f * phi_f_2

    return dpdz_F
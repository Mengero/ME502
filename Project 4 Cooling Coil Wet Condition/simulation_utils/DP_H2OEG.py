import numpy as np

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
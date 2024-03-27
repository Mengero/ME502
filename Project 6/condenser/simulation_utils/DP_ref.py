from condenser.simulation_utils.dpdz_f_2ph_Lopez_2014 import dpdz_f_2ph_Lopez_2014
from condenser.simulation_utils.dpdz_f_roundtube_churchill import dpdz_f_roundtube_churchill
from condenser.simulation_utils.void_fraction import void_fraction
import CoolProp.CoolProp as CP

def DP_ref_1ph(Re,roughness,rho_ref,G_ref,Vel_i,Vel_o,z_i,z_o,dx,Dh):
    '''
        for single phase pressure drop
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
    dpdz_f = dpdz_f_roundtube_churchill(Re,Dh,roughness)
    DP_static = dpdz_f * dx/Dh * G_ref**2/(2*rho_ref)

    # hydrostatic pressure: gravity, rho g h
    DP_hydrostatic = rho_ref*g*(z_i-z_o)
    # dynamic pressure: velocity, 1/2 rho v^2, for single phase, dynamic pressure is negligible
    # DP_dynamic = 1/2*rho_ref*(Vel_i**2-Vel_o**2)
    DP_dynamic = 0
    # total pressure=static + hydrostatic + dynamic
    DP_total = DP_static+DP_hydrostatic+DP_dynamic

    return DP_static,DP_hydrostatic,DP_dynamic,DP_total

def DP_ref_2ph(self,Re,roughness,P_ri,G_ref,x_ri,h_ro,z_i,z_o,dx,Dh):
    '''
        for two phase pressure drop
        * Input Variables:
            @Re             : Reynolds number
            @roughness      : roughness of inner tube [m]
            @P_ri           : inlet pressure, [kPa]
            @G_ref          : mass velocity of refriegrant [kg/m^2-s]
            @
            @z_i            : inlet height [m]
            @z_o            : outlet height [m]
            @dx             : length of each element [m]
            @Dh             : hydraulic diameter, [m]

         * Output Variables
            @DP_total       : [Pa]
    '''
    rho_ref=CP.PropsSI('D','P',P_ri*1e3,'Q',x_ri,self.fld)
    g=9.81          # [m/s**2]
    # static pressure: friction, (head loss)rho*g*h_L=rho*g*h_L,major + rho*g*h_L_minor
    #   rho*g*h_L_major=f*l/D*rho*V**2/2; rho*g*h_L_minor=K_L*rho*V**2/2
    #   negligible for minor one, only consider the major one
    T_sat = CP.PropsSI('T','P',P_ri*1e3,'Q',0.5,self.fld)-273.15
    dpdz_f = dpdz_f_2ph_Lopez_2014(self.fld,G_ref,T_sat,x_ri,Dh)        # [Pa/m]
    DP_static = dpdz_f * dx

    # hydrostatic pressure: gravity, rho g h
    DP_hydrostatic = rho_ref*g*(z_i-z_o)
    # dynamic pressure: velocity, 1/2 rho v^2, for single phase, dynamic pressure is negligible
    # DP_dynamic = 1/2*rho_ref*(Vel_i**2-Vel_o**2)
    P_ro_tmp=P_ri-DP_static
    h_ro_v_tmp=CP.PropsSI('H','P',P_ro_tmp*1e3,'Q',1,self.fld)/1e3      # [kJ/kg]
    h_ro_l_tmp=CP.PropsSI('H','P',P_ro_tmp*1e3,'Q',0,self.fld)/1e3      # [kJ/kg]
    x_ro_tmp=(h_ro-h_ro_l_tmp)/(h_ro_v_tmp-h_ro_l_tmp)
    rho_v_ri=CP.PropsSI('D','P',P_ri*1e3,'Q',1,self.fld)                # [kg/m^3]
    rho_l_ri=CP.PropsSI('D','P',P_ri*1e3,'Q',0,self.fld)
    rho_v_ro=CP.PropsSI('D','P',P_ro_tmp*1e3,'Q',1,self.fld)
    rho_l_ro=CP.PropsSI('D','P',P_ro_tmp*1e3,'Q',0,self.fld)
    [alpha_ri,rho_ref_tmp]=void_fraction(self,P_ri,x_ri)
    [alpha_ro_tmp,rho_ref_tmp]=void_fraction(self,P_ro_tmp,x_ro_tmp)
    A = x_ri**2/(rho_v_ri*alpha_ri)+(1-x_ri)**2/(rho_l_ri*(1-alpha_ri))
    B = x_ro_tmp**2/(rho_v_ro*alpha_ro_tmp)+(1-x_ro_tmp)**2/(rho_l_ro*(1-alpha_ro_tmp))
    DP_dynamic = G_ref**2*(B-A)
    # total pressure=static + hydrostatic + dynamic
    DP_total = DP_static+DP_hydrostatic+DP_dynamic

    return DP_static,DP_hydrostatic,DP_dynamic,DP_total
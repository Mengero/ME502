import numpy as np
import CoolProp.CoolProp as CP

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
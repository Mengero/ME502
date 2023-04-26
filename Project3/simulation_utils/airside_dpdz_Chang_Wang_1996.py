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
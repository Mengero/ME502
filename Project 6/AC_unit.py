import numpy as np
import CoolProp.CoolProp as CP
from compressor import compressor
from condenser.condenser import condenser,condenser_system
from evaporator.evaporator import evaporator,evaporator_system

def result_organizer(result_tmp):
    '''
    variables in result_tmp:
        'm_dot_ref','Cooling Capacity','COP','T_eao','Q_cond','Q_evap',
        'DP_r_cond','DP_r_evap','DP_a_cond','DP_a_evap',
        'SHR_evap','T_r_cond_dist','P_r_cond_dist','T_ao_cond_dist','T_ao_evap_dist',
        'Q_evap_dist','Q_cond_dist','Q_sen_evap_dist','Q_lat_evap_dist'
    '''
    result = {}
    variables = ['m_dot_ref','Cooling Capacity','COP','T_eao','Q_cond','Q_evap',\
                 'DP_r_cond','DP_r_evap','DP_a_cond','DP_a_evap',\
                 'SHR_evap','T_r_cond_dist','P_r_cond_dist','T_r_evap_dist','P_r_evap_dist',\
                 'T_ao_cond_dist','T_ao_evap_dist','Q_evap_dist','Q_cond_dist','Q_sen_evap_dist',\
                 'Q_lat_evap_dist']
    for i, value in enumerate(result_tmp):
            result[variables[i]] = (value)

    return result

def AC_unit_system(d_cp,eta_cp,eta_cpother,SC_cro,SH_ero,T_cai,P_cai,u_face_cai,T_eai,P_eai,phi_eai,m_dot_eai):
    '''
        This function integrate condenser and evaporator with pre-required SH and SC levels.
        Start with two guessing values of P_cpri and P_cpro, which corresponds to the evaporator outlet and 
        condenser inlet.

         * Input Variables:
            @d_cp           : compressor displacement, [m^3/hr]
            @eta_cp         : compressor efficiency, [-]
            @eta_cpother    : compressor other efficiency (motor and mechanical efficiency), [-]
            @SC_cro         : refriegrant subcool at the condenser outlet, [C]
            @SH_ero         : refriegrant superheat at the evaporator outlet, [C]
            @T_cai          : condenser air inlet temperature, [C]
            @P_cai          : condenser air inlet pressure, [kPa]
            @u_face_cai     : condenser air face velocity, [m/s]
            @T_eai          : evaporator air inlet temperature, [C]
            @P_eai          : evaporator air inlet pressure, [kPa]
            @phi_eai        : evaporator air inlet relative humidity, [-]
            @m_dot_eai      : evaporator air inlet mass flow rate, [kg/min]

         * Output Variables:
            @
    '''
    Ref = 'R1234yf'
    N_segment = 40              # [-], number of segment in each tube
    # Initial guessing of P_cpri and P_cpro
    P_cpro_tmp = 1200                   # [kPa]
    P_cpri_tmp = 500             # [kPa]

    DIFF_hxv = 1000          # [J/kg], enthalpy difference between expansion valve inlet and outlet
    DIFF_SC = 2             # [C], SC difference between predicted value and given value
    T_cpro_sat = CP.PropsSI('T','P',P_cpro_tmp*1e3,'Q',0.5,Ref)
    dPdT_cpro = (CP.PropsSI('P','T',T_cpro_sat+1,'Q',0.5,Ref)/1e3-P_cpro_tmp)/18
    dPdh_cpri = 2e-4

    iteration_count = 0

    while abs(DIFF_hxv)>50:
        iteration_count+=1
        ## ------------------------------- Compressor -------------------------------- ##
        T_cpri_sat = CP.PropsSI('T','P',P_cpri_tmp*1e3,'Q',0.5,Ref)-273.15          # [C]
        T_cpri = T_cpri_sat+SH_ero
        h_cpro,m_dot_ref = compressor(Ref,eta_cp,eta_cpother,P_cpri_tmp,P_cpro_tmp,T_cpri,d_cp)         # [J/kg], [kg/hr]
        m_dot_ref = m_dot_ref/3600      # [kg/hr] -> [kg/s]
        DIFF_SC = 5
        while abs(DIFF_SC)>0.2:
            ## ------------------------------- Compressor -------------------------------- ##
            T_cpri_sat = CP.PropsSI('T','P',P_cpri_tmp*1e3,'Q',0.5,Ref)-273.15          # [C]
            T_cpri = T_cpri_sat+SH_ero
            h_cpro,m_dot_ref = compressor(Ref,eta_cp,eta_cpother,P_cpri_tmp,P_cpro_tmp,T_cpri,d_cp)         # [J/kg], [kg/hr]
            m_dot_ref = m_dot_ref/3600      # [kg/hr] -> [kg/s] 
            ## ------------------------------- Condenser -------------------------------- ##
            h_cri = h_cpro; P_cri = P_cpro_tmp
            T_cri = CP.PropsSI('T','P',P_cri*1e3,'H',h_cri,Ref)-273.15      # [C]
            T_sat_cri = CP.PropsSI('T','P',P_cri*1e3,'Q',0.5,Ref)-273.15     # [C]
            htc_correlation = 'Shah'
            condenser_geoinfo = condenser(T_cri,T_sat_cri,m_dot_ref,T_cai,P_cai,u_face_cai,N_segment)
            condenser_result = condenser_system(condenser_geoinfo,htc_correlation)
            P_cro = condenser_result['P_ro']; T_cro = condenser_result['T_ro']
            T_sat_cro = CP.PropsSI('T','P',P_cro*1e3,'Q',0.5,Ref)-273.15
            if T_sat_cro<=T_cro:
                print('SC not efficient, recheck')
            h_cro = CP.PropsSI('H','P',P_cro*1e3,'T',T_cro+273.15,Ref)
            SC_tmp = T_sat_cro - T_cro          # positive value
            DIFF_SC = SC_tmp - SC_cro
            P_cpro_tmp -= DIFF_SC*dPdT_cpro       
        ## ------------------------------- Evaporator -------------------------------- ##
        m_dot_ref_tmp = m_dot_ref*1e3           # [g/s]
        evaporator_result = evaporator_system(Ref,T_cro,P_cro,m_dot_ref_tmp,SH_ero,T_eai,P_eai,phi_eai,m_dot_eai,N_segment,P_cpri_tmp)
        T_eri = evaporator_result['T_ri']; P_eri = evaporator_result['P_ri']; x_eri = evaporator_result['x_ri']
        h_eri = CP.PropsSI('H','P',P_eri*1e3,'Q',x_eri,Ref)
        DIFF_hxv = h_eri-h_cro
        P_cpri_tmp -= DIFF_hxv*dPdh_cpri
        print('------------------- Iteration # {} -------------------'.format(iteration_count))
        print('P_cpri = {} [kPa]'.format(P_cpri_tmp))
        print('P_cpro = {} [kPa]'.format(P_cpro_tmp))
        print('DIFF_hxv = {} [J/kg]'.format(DIFF_hxv))
        print('DIFF_SC = {} [C]'.format(DIFF_SC))
        print('m_dot_ref = {}'.format(m_dot_ref))

    print('iteartion Finished!')


    # ------------------------------- Process Output Result --------------------------------
    Cooling_Capacity = evaporator_result['Q_tot']/1e3; Q_cond = condenser_result['Q_tot']; Q_evap = Cooling_Capacity
    COP = Q_evap/(Q_cond-Q_evap); T_eao = evaporator_result['T_ao']
    DP_r_evap = evaporator_result['P_ri'] - evaporator_result['P_ro']; DP_a_evap = P_eai-evaporator_result['P_ao']
    DP_r_cond = P_cri - condenser_result['P_ro']; DP_a_cond = P_cai - condenser_result['P_ao']
    Q_sen_tot = np.sum(evaporator_result['Q_a_sen_dist'])*25/1e3
    SHR_evap = Q_sen_tot/Q_evap
    T_r_evap_dist = evaporator_result['T_r_dist']; P_r_evap_dist = evaporator_result['P_r_dist']
    T_ao_evap_dist = evaporator_result['T_a_dist']; T_edpa_dist = evaporator_result['Tdp_a_dist']
    T_r_cond_dist = condenser_result['T_r_dist']; P_r_cond_dist = condenser_result['P_r_dist']
    T_ao_cond_dist = condenser_result['T_a_dist']; Q_evap_dist = evaporator_result['Q_dist']
    Q_cond_dist = condenser_result['Q_dist']; Q_sen_evap_dist = evaporator_result['Q_a_sen_dist']
    Q_lat_evap_dist = evaporator_result['Q_a_lat_dist']
    

    result_tmp = [m_dot_ref,Cooling_Capacity,COP,T_eao,Q_cond,Q_evap,\
                  DP_r_cond,DP_r_evap,DP_a_cond,DP_a_evap,\
                  SHR_evap,T_r_cond_dist,P_r_cond_dist,T_r_evap_dist,P_r_evap_dist,\
                  T_ao_cond_dist,T_ao_evap_dist,Q_evap_dist,Q_cond_dist,\
                  Q_sen_evap_dist,Q_lat_evap_dist]
    result = result_organizer(result_tmp)

    return result


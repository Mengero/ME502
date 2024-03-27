from AC_unit import AC_unit_system
import pickle
import numpy as np

# ----------------------------------- Base Condition ----------------------------------- #
# d_cp = 5.5; eta_cp = 0.7; eta_cpother = 0.95; SC_cro = 5; SH_ero = 8; T_cai = 35; P_cai = 99.5; u_face_cai = 3
# T_eai = 35; P_eai = 99.5; phi_eai = 0.4; m_dot_eai = 9

# result_baseline = AC_unit_system(d_cp,eta_cp,eta_cpother,SC_cro,SH_ero,T_cai,P_cai,u_face_cai,T_eai,P_eai,phi_eai,m_dot_eai)
# with open("baseline_performance.pkl", "wb") as f:
#     pickle.dump((result_baseline), f)

# ----------------------------------- Effect of eta_cp Condition ----------------------------------- #
'''
    eta_cp: 0.6 ~ 0.8
'''
d_cp = 5.5; eta_cp = 0.7; eta_cpother = 0.95; SC_cro = 5; SH_ero = 8; T_cai = 35; P_cai = 99.5; u_face_cai = 3
T_eai = 35; P_eai = 99.5; phi_eai = 0.4; m_dot_eai = 9
result_etacp = {}

for eta_cp in np.linspace(0.6,0.8,num=10):
    print('--------------------------------------------------------------------------')
    print('solving cases for eta_cp = {}'.format(eta_cp))
    result_tmp = AC_unit_system(d_cp,eta_cp,eta_cpother,SC_cro,SH_ero,T_cai,P_cai,u_face_cai,T_eai,P_eai,phi_eai,m_dot_eai)
    result_etacp[eta_cp] = result_tmp

with open("eta_cp_impact.pkl", "wb") as f:
    pickle.dump((result_etacp), f)


# # ----------------------------------- Effect of SC_cro Condition ----------------------------------- #
# '''
#     SC_cro: 2C ~ 12C
# '''
# d_cp = 5.5; eta_cp = 0.7; eta_cpother = 0.95; SC_cro = 5; SH_ero = 8; T_cai = 35; P_cai = 99.5; u_face_cai = 3
# T_eai = 35; P_eai = 99.5; phi_eai = 0.4; m_dot_eai = 9
# result_SC_cro = {}

# for SC_cro in np.linspace(2,12,num=10):
#     print('--------------------------------------------------------------------------')
#     print('solving cases for SC_cro = {}'.format(SC_cro))
#     result_tmp = AC_unit_system(d_cp,eta_cp,eta_cpother,SC_cro,SH_ero,T_cai,P_cai,u_face_cai,T_eai,P_eai,phi_eai,m_dot_eai)
#     result_SC_cro[SC_cro] = result_tmp

# with open("SC_cro_impact.pkl", "wb") as f:
#     pickle.dump((result_SC_cro), f)


# ----------------------------------- Effect of SH_ero Condition ----------------------------------- #
# '''
#     SH_ero: 2C ~ 12C
# '''
# d_cp = 5.5; eta_cp = 0.7; eta_cpother = 0.95; SC_cro = 5; SH_ero = 8; T_cai = 35; P_cai = 99.5; u_face_cai = 3
# T_eai = 35; P_eai = 99.5; phi_eai = 0.4; m_dot_eai = 9
# result_SH_ero = {}

# for SH_ero in np.linspace(2,12,num=10):
#     print('--------------------------------------------------------------------------')
#     print('solving cases for SH_ero = {}'.format(SH_ero))
#     result_tmp = AC_unit_system(d_cp,eta_cp,eta_cpother,SC_cro,SH_ero,T_cai,P_cai,u_face_cai,T_eai,P_eai,phi_eai,m_dot_eai)
#     result_SH_ero[SH_ero] = result_tmp

# with open("SH_ero_impact.pkl", "wb") as f:
#     pickle.dump((result_SH_ero), f)


# ----------------------------------- Effect of T_cai Condition ----------------------------------- #
# '''
#     T_cai: 30C ~ 40C
# '''
# d_cp = 5.5; eta_cp = 0.7; eta_cpother = 0.95; SC_cro = 5; SH_ero = 8; T_cai = 35; P_cai = 99.5; u_face_cai = 3
# T_eai = 35; P_eai = 99.5; phi_eai = 0.4; m_dot_eai = 9
# result_T_cai = {}

# for T_cai in np.linspace(30,40,num=10):
#     print('--------------------------------------------------------------------------')
#     print('solving cases for T_cai = {}'.format(T_cai))
#     result_tmp = AC_unit_system(d_cp,eta_cp,eta_cpother,SC_cro,SH_ero,T_cai,P_cai,u_face_cai,T_eai,P_eai,phi_eai,m_dot_eai)
#     result_T_cai[T_cai] = result_tmp

# with open("T_cai_impact.pkl", "wb") as f:
#     pickle.dump((result_T_cai), f)

# # ----------------------------------- Effect of T_eai Condition ----------------------------------- #
# '''
#     T_eai: 30C ~ 40C
# '''
# d_cp = 5.5; eta_cp = 0.7; eta_cpother = 0.95; SC_cro = 5; SH_ero = 8; T_cai = 35; P_cai = 99.5; u_face_cai = 3
# T_eai = 35; P_eai = 99.5; phi_eai = 0.4; m_dot_eai = 9
# result_T_eai = {}

# for T_eai in np.linspace(30,40,num=10):
#     print('--------------------------------------------------------------------------')
#     print('solving cases for T_eai = {}'.format(T_eai))
#     result_tmp = AC_unit_system(d_cp,eta_cp,eta_cpother,SC_cro,SH_ero,T_cai,P_cai,u_face_cai,T_eai,P_eai,phi_eai,m_dot_eai)
#     result_T_eai[T_eai] = result_tmp

# with open("T_eai_impact.pkl", "wb") as f:
#     pickle.dump((result_T_eai), f)

# ----------------------------------- Effect of phi_eai Condition ----------------------------------- #
'''
    phi_eai: 0.3 ~ 0.5
'''
d_cp = 5.5; eta_cp = 0.7; eta_cpother = 0.95; SC_cro = 5; SH_ero = 8; T_cai = 35; P_cai = 99.5; u_face_cai = 3
T_eai = 35; P_eai = 99.5; phi_eai = 0.4; m_dot_eai = 9
result_phi_eai = {}

for phi_eai in np.linspace(0.3,0.5,num=10):
    print('--------------------------------------------------------------------------')
    print('solving cases for phi_eai = {}'.format(phi_eai))
    result_tmp = AC_unit_system(d_cp,eta_cp,eta_cpother,SC_cro,SH_ero,T_cai,P_cai,u_face_cai,T_eai,P_eai,phi_eai,m_dot_eai)
    result_phi_eai[phi_eai] = result_tmp

with open("phi_eai_impact.pkl", "wb") as f:
    pickle.dump((result_phi_eai), f)
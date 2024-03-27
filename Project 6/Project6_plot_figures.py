import pickle
import matplotlib.pyplot as plt
import numpy as np

'''
m_dot_ref,Cooling_Capacity,COP,T_eao,Q_cond,Q_evap,\
DP_r_cond_dist,DP_r_evap_dist,DP_a_cond_dist,DP_a_evap_dist,\
SHR_evap,T_r_cond_dist,P_r_cond_dist,T_r_evap_dist,P_r_evap_dist,\
T_ao_cond_dist,T_ao_evap_dist,Q_evap_dist,Q_cond_dist,\
Q_sen_evap_dist,Q_lat_evap_dist
'''


N_elemptube = 161
with open("baseline_performance.pkl", "rb") as f:
    result_baseline = pickle.load(f)
with open("eta_cp_impact.pkl", "rb") as f:
    result_eta_cp = pickle.load(f)
with open("SC_cro_impact.pkl", "rb") as f:
    result_SC_cro = pickle.load(f)
with open("SH_ero_impact.pkl", "rb") as f:
    result_SH_ero = pickle.load(f)
with open("T_cai_impact.pkl", "rb") as f:
    result_T_cai = pickle.load(f)
with open("T_eai_impact.pkl", "rb") as f:
    result_T_eai = pickle.load(f)
with open("phi_eai_impact.pkl", "rb") as f:
    result_phi_eai = pickle.load(f)
folder_path = 'Figures/'
# --------------------------------------------------------- Overall Performance --------------------------------------------------------- #
x = range(1,N_elemptube)


fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['T_r_cond_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig(folder_path+'Overall_T_r_cond.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['T_r_evap_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig(folder_path+'Overall_T_r_evap.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['P_r_cond_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('P [kPa]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig(folder_path+'Overall_P_r_cond.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['P_r_evap_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('P [kPa]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig(folder_path+'Overall_P_r_evap.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['T_ao_cond_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig(folder_path+'Overall_T_ao_cond.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['T_ao_evap_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig(folder_path+'Overall_T_ao_evap.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(x,result_baseline['Q_evap_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('Q [W]',fontsize=20)
ax1.set_xlabel('Circuit Location',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'Overall_Q_evap.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(x,result_baseline['Q_cond_dist']*1e3)
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('Q [W]',fontsize=20)
ax1.set_xlabel('Circuit Location',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'Overall_Q_cond.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(x,result_baseline['Q_sen_evap_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('Q [W]',fontsize=20)
ax1.set_xlabel('Circuit Location',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'Overall_Q_sen_evap.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(x,result_baseline['Q_lat_evap_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
ax1.set_ylabel('Q [W]',fontsize=20)
ax1.set_xlabel('Circuit Location',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'Overall_Q_lat_evap.png')


# --------------------------------------------------------- Impact of eta_cp --------------------------------------------------------- #
m_dot_ref_tmp=[];CC_tmp=[];COP_tmp=[];T_eao=[]
for key in np.linspace(0.6,0.8,num=10):
    m_dot_ref_tmp.append(result_eta_cp[key]['m_dot_ref']*1e3)
    CC_tmp.append(result_eta_cp[key]['Cooling Capacity'])
    COP_tmp.append(result_eta_cp[key]['COP'])
    T_eao.append(result_eta_cp[key]['T_eao'])
eta_cp = np.linspace(0.6,0.8,num=10)

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(eta_cp,m_dot_ref_tmp)
ax1.grid()
ax1.set_ylabel('mass flow rate [g/s]',fontsize=20)
ax1.set_xlabel('eta_cp [-]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'eta_cp_MFR.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(eta_cp,CC_tmp)
ax1.grid()
ax1.set_ylabel('Cooling Capacity [kW]',fontsize=20)
ax1.set_xlabel('eta_cp [-]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'eta_cp_CC.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(eta_cp,COP_tmp)
ax1.grid()
ax1.set_ylabel('COP [-]',fontsize=20)
ax1.set_xlabel('eta_cp [-]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'eta_cp_COP.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(eta_cp,T_eao)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=20)
ax1.set_xlabel('eta_cp [-]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'eta_cp_Teao.png')

# --------------------------------------------------------- Impact of SC_cro --------------------------------------------------------- #
m_dot_ref_tmp=[];CC_tmp=[];COP_tmp=[];T_eao=[]
for key in np.linspace(2,12,num=10):
    m_dot_ref_tmp.append(result_SC_cro[key]['m_dot_ref']*1e3)
    CC_tmp.append(result_SC_cro[key]['Cooling Capacity'])
    COP_tmp.append(result_SC_cro[key]['COP'])
    T_eao.append(result_SC_cro[key]['T_eao'])
SC_cro = np.linspace(2,12,num=10)

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(SC_cro,m_dot_ref_tmp)
ax1.grid()
ax1.set_ylabel('mass flow rate [g/s]',fontsize=20)
ax1.set_xlabel('SC_cro [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'SC_cro_MFR.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(SC_cro,CC_tmp)
ax1.grid()
ax1.set_ylabel('Cooling Capacity [kW]',fontsize=20)
ax1.set_xlabel('SC_cro [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'SC_cro_CC.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(SC_cro,COP_tmp)
ax1.grid()
ax1.set_ylabel('COP [-]',fontsize=20)
ax1.set_xlabel('SC_cro [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'SC_cro_COP.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(SC_cro,T_eao)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=20)
ax1.set_xlabel('SC_cro [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'SC_cro_Teao.png')

# --------------------------------------------------------- Impact of SH_ero --------------------------------------------------------- #
m_dot_ref_tmp=[];CC_tmp=[];COP_tmp=[];T_eao=[]
for key in np.linspace(2,12,num=10):
    m_dot_ref_tmp.append(result_SH_ero[key]['m_dot_ref']*1e3)
    CC_tmp.append(result_SH_ero[key]['Cooling Capacity'])
    COP_tmp.append(result_SH_ero[key]['COP'])
    T_eao.append(result_SH_ero[key]['T_eao'])
SH_ero = np.linspace(2,12,num=10)

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(SH_ero,m_dot_ref_tmp)
ax1.grid()
ax1.set_ylabel('mass flow rate [g/s]',fontsize=20)
ax1.set_xlabel('SH_ero [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'SH_ero_MFR.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(SH_ero,CC_tmp)
ax1.grid()
ax1.set_ylabel('Cooling Capacity [kW]',fontsize=20)
ax1.set_xlabel('SH_ero [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'SH_ero_CC.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(SH_ero,COP_tmp)
ax1.grid()
ax1.set_ylabel('COP [-]',fontsize=20)
ax1.set_xlabel('SH_ero [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'SH_ero_COP.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(SH_ero,T_eao)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=20)
ax1.set_xlabel('SH_ero [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'SH_ero_Teao.png')

# --------------------------------------------------------- Impact of T_cai --------------------------------------------------------- #
m_dot_ref_tmp=[];CC_tmp=[];COP_tmp=[];T_eao=[]
for key in np.linspace(30,40,num=10):
    m_dot_ref_tmp.append(result_T_cai[key]['m_dot_ref']*1e3)
    CC_tmp.append(result_T_cai[key]['Cooling Capacity'])
    COP_tmp.append(result_T_cai[key]['COP'])
    T_eao.append(result_T_cai[key]['T_eao'])
T_cai = np.linspace(30,40,num=10)

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(T_cai,m_dot_ref_tmp)
ax1.grid()
ax1.set_ylabel('mass flow rate [g/s]',fontsize=20)
ax1.set_xlabel('T_cai [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'T_cai_MFR.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(T_cai,CC_tmp)
ax1.grid()
ax1.set_ylabel('Cooling Capacity [kW]',fontsize=20)
ax1.set_xlabel('T_cai [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'T_cai_CC.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(T_cai,COP_tmp)
ax1.grid()
ax1.set_ylabel('COP [-]',fontsize=20)
ax1.set_xlabel('T_cai [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'T_cai_COP.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(T_cai,T_eao)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=20)
ax1.set_xlabel('T_cai [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'T_cai_Teao.png')

# --------------------------------------------------------- Impact of T_eai --------------------------------------------------------- #
m_dot_ref_tmp=[];CC_tmp=[];COP_tmp=[];T_eao=[]
for key in np.linspace(30,40,num=10):
    m_dot_ref_tmp.append(result_T_eai[key]['m_dot_ref']*1e3)
    CC_tmp.append(result_T_eai[key]['Cooling Capacity'])
    COP_tmp.append(result_T_eai[key]['COP'])
    T_eao.append(result_T_eai[key]['T_eao'])
T_eai = np.linspace(30,40,num=10)

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(T_eai,m_dot_ref_tmp)
ax1.grid()
ax1.set_ylabel('mass flow rate [g/s]',fontsize=20)
ax1.set_xlabel('T_eai [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'T_eai_MFR.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(T_eai,CC_tmp)
ax1.grid()
ax1.set_ylabel('Cooling Capacity [kW]',fontsize=20)
ax1.set_xlabel('T_eai [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'T_eai_CC.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(T_eai,COP_tmp)
ax1.grid()
ax1.set_ylabel('COP [-]',fontsize=20)
ax1.set_xlabel('T_eai [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'T_eai_COP.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(T_eai,T_eao)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=20)
ax1.set_xlabel('T_eai [C]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'T_eai_Teao.png')

# --------------------------------------------------------- Impact of phi_eai --------------------------------------------------------- #
m_dot_ref_tmp=[];CC_tmp=[];COP_tmp=[];T_eao=[]
for key in np.linspace(0.3,0.5,num=10):
    m_dot_ref_tmp.append(result_phi_eai[key]['m_dot_ref']*1e3)
    CC_tmp.append(result_phi_eai[key]['Cooling Capacity'])
    COP_tmp.append(result_phi_eai[key]['COP'])
    T_eao.append(result_phi_eai[key]['T_eao'])
phi_eai = np.linspace(0.3,0.5,num=10)

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(phi_eai,m_dot_ref_tmp)
ax1.grid()
ax1.set_ylabel('mass flow rate [g/s]',fontsize=20)
ax1.set_xlabel('phi_eai [%]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'phi_eai_MFR.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(phi_eai,CC_tmp)
ax1.grid()
ax1.set_ylabel('Cooling Capacity [kW]',fontsize=20)
ax1.set_xlabel('phi_eai [%]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'phi_eai_CC.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(phi_eai,COP_tmp)
ax1.grid()
ax1.set_ylabel('COP [-]',fontsize=20)
ax1.set_xlabel('phi_eai [%]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'phi_eai_COP.png')

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(phi_eai,T_eao)
ax1.grid()
ax1.set_ylabel('T [C]',fontsize=20)
ax1.set_xlabel('phi_eai [%]',fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)
fig.savefig(folder_path+'phi_eai_Teao.png')
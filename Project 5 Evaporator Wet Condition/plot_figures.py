import pickle
import matplotlib.pyplot as plt
import numpy as np
N_elemptube = 161
with open("baseline_performance.pkl", "rb") as f:
    result_baseline = pickle.load(f)
with open("RH_ai_impact.pkl", "rb") as f:
    result_RH_ai = pickle.load(f)
# ----------------------------- Overall Performance ---------------------------- #
Q_sen_tot = np.sum(result_baseline['Q_a_sen_dist'])*25/1e3
Q_lat_tot = np.sum(result_baseline['Q_a_lat_dist'])*25/1e3
Q_sen_tot_02=np.sum(result_RH_ai[0.2]['Q_a_sen_dist'])*25/1e3
Q_lat_tot_02=np.sum(result_RH_ai[0.2]['Q_a_lat_dist'])*25/1e3
Q_sen_tot_08=np.sum(result_RH_ai[0.8]['Q_a_sen_dist'])*25/1e3
Q_lat_tot_08=np.sum(result_RH_ai[0.8]['Q_a_lat_dist'])*25/1e3
x = range(1,N_elemptube)
fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['T_r_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('EG Temperature Plot Through the Tube, N_segment=40')
ax1.set_ylabel('T [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('Overall_T_r.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['T_ai_dist'],label='air inlet temperature')
ax1.plot(x,result_baseline['T_a_dist'],label='air outlet temperature')
ax1.set_xlim(-1,N_elemptube)
ax1.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)
ax1.grid()
# ax1.set_title('Air Temperature Plot Through the Tube, N_segment=40')
ax1.set_ylabel('T [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('Overall_T_a.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['Tdp_a_dist'],label='air outlet Tdp')
ax1.set_xlim(-1,N_elemptube)
ax1.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)
ax1.grid()
# ax1.set_title('Air Dew Point Temperature Plot Through the Tube, N_segment=40')
ax1.set_ylabel('Tdp [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('Overall_Tdp_a.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_baseline['P_r_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('EG Pressure Plot Through the Tube, N_segment=40')
ax1.set_ylabel('P [kPa]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('Overall_P_r.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(x,result_baseline['htc_r_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('EG HTC Plot Through the Tube, N_segment=40')
ax1.set_ylabel('HTC [W/m^2-K]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('Overall_HTC_r.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(x,result_baseline['Q_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Heat Transfer Rate Plot Through the Tube, N_segment=40')
ax1.set_ylabel('Q [W]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('Overall_Q.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(x,result_baseline['Q_a_sen_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Sensible Heat Transfer Rate Plot Through the Tube, N_segment=40')
ax1.set_ylabel('Q [W]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('Overall_Q_sen.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(x,result_baseline['Q_a_lat_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Latent Heat Transfer Rate Plot Through the Tube, N_segment=40')
ax1.set_ylabel('Q [W]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('Overall_Q_lat.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(x,result_baseline['SHR_dist'])
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Latent Heat Transfer Rate Plot Through the Tube, N_segment=40')
ax1.set_ylabel('SHR [%]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('Overall_SHR.png')

# ----------------------- Impact of RH_ai ------------------------ #
fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai]['T_r_dist'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Impact of RH_ai on EG Temperature, N_segment=40')
ax1.set_ylabel('T [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('RH_ai_T_r.png')

fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai]['T_a_dist'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Impact of RH_ai on Air Outlet Temperature, N_segment=40')
ax1.set_ylabel('T [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('RH_ai_T_ao.png')

fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai]['Tdp_a_dist'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Impact of RH_ai on Air Outlet Dew Point Temperature, N_segment=40')
ax1.set_ylabel('Tdp [C]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('RH_ai_Tdp_ao.png')

fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai]['P_r_dist'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Impact of RH_ai on EG Pressure, N_segment=40')
ax1.set_ylabel('P [kPa]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('RH_ai_P_r.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai]['htc_r_dist'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Impact of RH_ai on EG HTC, N_segment=40')
ax1.set_ylabel('HTC [W/m^2-K]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('RH_ai_htc_r.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai]['Q_dist'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Impact of RH_ai on segment heat transfer rate, N_segment=40')
ax1.set_ylabel('Q [W]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('RH_ai_Q.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai]['Q_a_sen_dist'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Impact of RH_ai on segment sensible heat transfer rate, N_segment=40')
ax1.set_ylabel('Q [W]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('RH_ai_Q_sen.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai]['Q_a_lat_dist'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Impact of RH_ai on segment latent heat transfer rate, N_segment=40')
ax1.set_ylabel('Q [W]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('RH_ai_Q_latent.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai]['SHR_dist'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,N_elemptube)
ax1.grid()
# ax1.set_title('Impact of RH_ai on segment latent heat transfer rate, N_segment=40')
ax1.set_ylabel('SHR [%]',fontsize=16)
ax1.set_xlabel('Circuit Location',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)
fig.savefig('RH_ai_SHR.png')

# ----------------- Generate Legend for RH_ai and vol_dot_a ----------------- #
x = []; y = []
fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,y,label='RH_ai = {}'.format(RH_ai))
plt.legend(loc='center',fontsize=12)
fig.savefig('RH_ai_legend.png')
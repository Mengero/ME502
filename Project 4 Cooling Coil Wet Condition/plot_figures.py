import pickle
import matplotlib.pyplot as plt
import numpy as np

with open("overall_performance.pkl", "rb") as f:
    result_overall = pickle.load(f)
with open("RH_ai_impact.pkl", "rb") as f:
    result_RH_ai = pickle.load(f)
with open("vol_dot_a_impact.pkl", "rb") as f:
    result_vol_dot_a = pickle.load(f)
# ----------------------------- Overall Performance ---------------------------- #
x = range(1,41)
fig, ax1 = plt.subplots()
ax1.plot(x,result_overall[0]['EG Temperature'])
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('EG Temperature Plot Through the Tube, N_segment=40')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Overall_T_r.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_overall[0]['air inlet temperature'],label='air inlet temperature')
ax1.plot(x,result_overall[0]['air outlet temperature'],label='air outlet temperature')
ax1.set_xlim(-1,41)
ax1.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)
ax1.grid()
# ax1.set_title('Air Temperature Plot Through the Tube, N_segment=40')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Overall_T_a.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_overall[0]['air inlet dew point temperature'],label='air inlet Tdp')
ax1.plot(x,result_overall[0]['air outlet dew point temperature'],label='air outlet Tdp')
ax1.set_xlim(-1,41)
ax1.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)
ax1.grid()
# ax1.set_title('Air Dew Point Temperature Plot Through the Tube, N_segment=40')
ax1.set_ylabel('Tdp [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Overall_Tdp_a.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_overall[0]['EG Pressure'])
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('EG Pressure Plot Through the Tube, N_segment=40')
ax1.set_ylabel('P [kPa]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Overall_P_r.png')

fig, ax1 = plt.subplots()
ax1.plot(x,result_overall[0]['EG HTC'])
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('EG HTC Plot Through the Tube, N_segment=40')
ax1.set_ylabel('HTC [kW/m^2-K]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Overall_HTC_r.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(x,result_overall[0]['segment heat transfer rate'])
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Heat Transfer Rate Plot Through the Tube, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Overall_Q.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(x,result_overall[0]['segment sensible heat transfer rate'])
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Sensible Heat Transfer Rate Plot Through the Tube, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Overall_Q_sen.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(x,result_overall[0]['segment latent heat transfer rate'])
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Latent Heat Transfer Rate Plot Through the Tube, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Overall_Q_lat.png')

# ----------------------- Impact of RH_ai ------------------------ #
fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai][0]['EG Temperature'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of RH_ai on EG Temperature, N_segment=40')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('RH_ai_T_r.png')

fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai][0]['air outlet temperature'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of RH_ai on Air Outlet Temperature, N_segment=40')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('RH_ai_T_ao.png')

fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai][0]['air outlet dew point temperature'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of RH_ai on Air Outlet Dew Point Temperature, N_segment=40')
ax1.set_ylabel('Tdp [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('RH_ai_Tdp_ao.png')

fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai][0]['EG Pressure'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of RH_ai on EG Pressure, N_segment=40')
ax1.set_ylabel('P [kPa]')
ax1.set_xlabel('Circuit Location')
fig.savefig('RH_ai_P_r.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai][0]['EG HTC'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of RH_ai on EG HTC, N_segment=40')
ax1.set_ylabel('HTC [kW/m^2-K]')
ax1.set_xlabel('Circuit Location')
fig.savefig('RH_ai_htc_r.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai][0]['segment heat transfer rate'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of RH_ai on segment heat transfer rate, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('RH_ai_Q.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai][0]['segment sensible heat transfer rate'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of RH_ai on segment sensible heat transfer rate, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('RH_ai_Q_sen.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,result_RH_ai[RH_ai][0]['segment latent heat transfer rate'],label='RH_ai={}'.format(RH_ai))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of RH_ai on segment latent heat transfer rate, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('RH_ai_Q_latent.png')

# ----------------------- Impact of vol_dot_a ------------------------ #
fig, ax1 = plt.subplots()
for vol_dot_a in np.linspace(250,600,num=10):
    ax1.plot(x,result_vol_dot_a[vol_dot_a][0]['EG Temperature'],label='vol_dot_a={}'.format(vol_dot_a))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of vol_dot_a on EG Temperature, N_segment=40')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('vol_dot_a_T_r.png')

fig, ax1 = plt.subplots()
for vol_dot_a in np.linspace(250,600,num=10):
    ax1.plot(x,result_vol_dot_a[vol_dot_a][0]['air outlet temperature'],label='vol_dot_a={}'.format(vol_dot_a))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of vol_dot_a on Air Outlet Temperature, N_segment=40')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('vol_dot_a_T_ao.png')

fig, ax1 = plt.subplots()
for vol_dot_a in np.linspace(250,600,num=10):
    ax1.plot(x,result_vol_dot_a[vol_dot_a][0]['air outlet dew point temperature'],label='vol_dot_a={}'.format(vol_dot_a))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of vol_dot_a on Air Outlet Dew Point Temperature, N_segment=40')
ax1.set_ylabel('Tdp [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('vol_dot_a_Tdp_ao.png')

fig, ax1 = plt.subplots()
for vol_dot_a in np.linspace(250,600,num=10):
    ax1.plot(x,result_vol_dot_a[vol_dot_a][0]['EG Pressure'],label='vol_dot_a={}'.format(vol_dot_a))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of vol_dot_a on EG Pressure, N_segment=40')
ax1.set_ylabel('P [kPa]')
ax1.set_xlabel('Circuit Location')
fig.savefig('vol_dot_a_P_r.png')

fig, ax1 = plt.subplots()
for vol_dot_a in np.linspace(250,600,num=10):
    ax1.plot(x,result_vol_dot_a[vol_dot_a][0]['EG HTC'],label='vol_dot_a={}'.format(vol_dot_a))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of vol_dot_a on EG HTC, N_segment=40')
ax1.set_ylabel('HTC [kW/m^2-K]')
ax1.set_xlabel('Circuit Location')
fig.savefig('vol_dot_a_htc_r.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for vol_dot_a in np.linspace(250,600,num=10):
    ax1.plot(x,result_vol_dot_a[vol_dot_a][0]['segment heat transfer rate'],label='vol_dot_a={}'.format(vol_dot_a))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of vol_dot_a on segment heat transfer rate, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('vol_dot_a_Q.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for vol_dot_a in np.linspace(250,600,num=10):
    ax1.plot(x,result_vol_dot_a[vol_dot_a][0]['segment sensible heat transfer rate'],label='vol_dot_a={}'.format(vol_dot_a))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of vol_dot_a on segment sensible heat transfer rate, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('vol_dot_a_Q_sen.png')

fig, ax1 = plt.subplots(figsize=(8, 6))
for vol_dot_a in np.linspace(250,600,num=10):
    ax1.plot(x,result_vol_dot_a[vol_dot_a][0]['segment latent heat transfer rate'],label='vol_dot_a={}'.format(vol_dot_a))
ax1.set_xlim(-1,41)
ax1.grid()
# ax1.set_title('Impact of vol_dot_a on segment latent heat transfer rate, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('vol_dot_a_Q_latent.png')

# ----------------- Generate Legend for RH_ai and vol_dot_a ----------------- #
x = []; y = []
fig, ax1 = plt.subplots()
for RH_ai in np.linspace(0.2,0.8,num=10):
    ax1.plot(x,y,label='RH_ai = {}'.format(RH_ai))
plt.legend(loc='center')
fig.savefig('RH_ai_legend.png')

fig, ax1 = plt.subplots()
for vol_dot_a in np.linspace(250,600,num=10):
    ax1.plot(x,y,label='vol_dot_a = {}'.format(vol_dot_a))
plt.legend(loc='center')
fig.savefig('vol_dot_a_legend.png')
import condenser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

'''
    This document finish the requirement in Project #2:
        1. Compare two different correlations for 2ph htc calculation
        2. Compare two different correlations for 2ph dpdz calculation
        3. Impact of element number on simulation results

'''

# ----------------- Condenser Input Variables --------------- #
fld = 'R1234yf'
T_ri = 75       # [C]
T_sat_ri = 48   # [C]
m_dot_cr = 35/1e3   # [g/s]
T_ai = 35       # [C]
P_ai = 99.5     # [kPa]
vol_dot_cai = 1500    # [CFM], cubic feet per minute
htc_correlation = 'Shah'
N_segment = 10
self=condenser.condenser(T_ri,T_sat_ri,m_dot_cr,T_ai,P_ai,vol_dot_cai,N_segment)


[T_ro,P_ro,x_ro,T_ao,P_ao,Q_tot,T_ro_dict,P_ro_dict,x_ro_dict,Q_dist_dict,T_r_dist_dict,T_a_dist_dict,P_r_dist_dict,P_a_dist_dict,DELTAP_air_dist_dict,DELTAP_ref_dist_dict,x_r_dist_dict,htc_ref_dist_dict] = condenser.condenser_system(self,htc_correlation)

# ---------------- Different number of segments impact on simulation results -------------- #
# Define data storage
N_segment_1_dict = {}
N_segment_2_dict = {}
Tro_dict = {}
Pro_dict = {}
Tao_dict = {}
Pao_dict = {}
xro_dict = {}
Qdist_dict = {}
Q_dict = {}
Tr_dist_dict = {}
Ta_dist_dict = {}
Pr_dist_dict = {}
Pa_dist_dict = {}
xr_dist_dict = {}
DELTAPair_dist_dict = {}
DELTAPref_dist_dict = {}
htcref_dist_dict = {}

for N_segment in range(1,80):
    N_segment_1 = np.array(range(N_segment+1))
    N_segment_2 = np.array(range(1,N_segment+1))
    self=condenser.condenser(T_ri,T_sat_ri,m_dot_cr,T_ai,P_ai,vol_dot_cai,N_segment)
    [T_ro,P_ro,x_ro,T_ao,P_ao,Q_tot,T_ro_dict,P_ro_dict,x_ro_dict,Q_dist_dict,T_r_dist_dict,T_a_dist_dict,P_r_dist_dict,P_a_dist_dict,DELTAP_air_dist_dict,DELTAP_ref_dist_dict,x_r_dist_dict,htc_ref_dist_dict] = condenser.condenser_system(self,htc_correlation)
    Q_dist = np.concatenate((Q_dist_dict[0],Q_dist_dict[1]),axis=None); T_r_dist = np.concatenate((T_r_dist_dict[0],T_r_dist_dict[1]),axis=None); P_r_dist = np.concatenate((P_r_dist_dict[0],P_r_dist_dict[1]),axis=None)
    T_a_dist = np.concatenate((T_a_dist_dict[0],T_a_dist_dict[1]),axis=None); P_a_dist = np.concatenate((P_a_dist_dict[0],P_a_dist_dict[1]),axis=None); x_r_dist = np.concatenate((x_r_dist_dict[0],x_r_dist_dict[1]),axis=None)
    DELTAP_air_dist = np.concatenate((DELTAP_air_dist_dict[0],DELTAP_air_dist_dict[1]),axis=None); DELTAP_ref_dist = np.concatenate((DELTAP_ref_dist_dict[0],DELTAP_ref_dist_dict[1]),axis=None)
    htc_ref_dist = np.concatenate((htc_ref_dist_dict[0],htc_ref_dist_dict[1]),axis=None)
    N_segment_1_dict[N_segment] = N_segment_1
    N_segment_2_dict[N_segment] = N_segment_2
    Tro_dict[N_segment] = T_ro; Pro_dict[N_segment] = P_ro
    Tao_dict[N_segment] = T_ao; Pao_dict[N_segment] = P_ao
    xro_dict[N_segment] = x_ro
    Qdist_dict[N_segment] = Q_dist
    Q_dict[N_segment] = Q_tot
    Tr_dist_dict[N_segment] = T_r_dist; Ta_dist_dict[N_segment] = T_a_dist
    Pr_dist_dict[N_segment] = P_r_dist; Pa_dist_dict[N_segment] = P_a_dist
    xr_dist_dict[N_segment] = x_r_dist
    DELTAPair_dist_dict[N_segment] = DELTAP_air_dist
    DELTAPref_dist_dict[N_segment] = DELTAP_ref_dist
    htcref_dist_dict[N_segment] = htc_ref_dist

# Plot Figures

# Tro
fig, ax1 = plt.subplots()
T_ro_list=Tro_dict.items()
x,y = zip(*T_ro_list)
ax1.plot(x,y)
ax1.grid()
ax1.set_title('Condenser Temperature Outlet')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Number of Segment [-]')
fig.savefig('Nsegment_Tro.pdf')
# Pro
fig, ax1 = plt.subplots()
P_ro_list=Pro_dict.items()
x,y = zip(*P_ro_list)
ax1.plot(x,y)
ax1.grid()
ax1.set_title('Condenser Pressure Outlet')
ax1.set_ylabel('P [kPa]')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Number of Segment [-]')
fig.savefig('Nsegment_Pro.pdf')
#xro
fig, ax1 = plt.subplots()
x_ro_list=xro_dict.items()
x,y = zip(*x_ro_list)
ax1.plot(x,y)
ax1.grid()
ax1.set_title('Condenser Vapor Quality Outlet')
ax1.set_ylabel('Vapor Quality [-]')
# ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Number of Segment [-]')
fig.savefig('Nsegment_xro.pdf')
#Tao
fig, ax1 = plt.subplots()
T_ao_list=Tao_dict.items()
x,y = zip(*T_ao_list)
ax1.plot(x,y)
ax1.grid()
ax1.set_title('Condenser Air Outlet Temperature')
ax1.set_ylabel('T [C]')
# ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Number of Segment [-]')
fig.savefig('Nsegment_Tao.pdf')
#Pao
fig, ax1 = plt.subplots()
P_ao_list=Pao_dict.items()
x,y = zip(*P_ao_list)
ax1.plot(x,y)
ax1.grid()
ax1.set_title('Condenser Air Outlet Pressure')
ax1.set_ylabel('P [kPa]')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Number of Segment [-]')
fig.savefig('Nsegment_Pao.pdf')
#Qtot
fig, ax1 = plt.subplots()
Q_list=Q_dict.items()
x,y = zip(*Q_list)
ax1.plot(x,y)
ax1.grid()
ax1.set_title('Condenser Capacity')
ax1.set_ylabel('Q [kW]')
# ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Number of Segment [-]')
fig.savefig('Nsegment_Qtot.pdf')

# Condenser Performance
#Tr
N_segment = 60
x=np.linspace(1,2*N_segment,num=2*N_segment)
fig, ax1 = plt.subplots()
ax1.plot(x,Tr_dist_dict[N_segment])
ax1.axvline(x=60,color='black',linestyle="--")
ax1.grid()
ax1.set_title('Condenser Temperature Distribution, N=60')
ax1.set_ylabel('T [C]')
# ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Loaction along tube (sgement) [-]')
fig.savefig('Performance_Tr.pdf')
#Pr
fig, ax1 = plt.subplots()
ax1.plot(x,Pr_dist_dict[N_segment])
ax1.axvline(x=60,color='black',linestyle="--")
ax1.grid()
ax1.set_title('Condenser Pressure Distribution, N=60')
ax1.set_ylabel('P [kPa]')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Loaction along tube (sgement) [-]')
fig.savefig('Performance_Pr.pdf')
#xr
fig, ax1 = plt.subplots()
ax1.plot(x,xr_dist_dict[N_segment])
ax1.axvline(x=60,color='black',linestyle="--")
ax1.grid()
ax1.set_title('Condenser Vapor Quality Distribution, N=60')
ax1.set_ylabel('Vapor Quality [-]')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Loaction along tube (sgement) [-]')
fig.savefig('Performance_xr.pdf')
#Tao
fig, ax1 = plt.subplots()
ax1.plot(x,Ta_dist_dict[N_segment])
ax1.axvline(x=60,color='black',linestyle="--")
ax1.grid()
ax1.set_title('Condenser Air Outlet Temperature Distribution, N=60')
ax1.set_ylabel('T [C]')
# ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Loaction along tube (sgement) [-]')
fig.savefig('Performance_Tao.pdf')
#Pao
fig, ax1 = plt.subplots()
ax1.plot(x,Pa_dist_dict[N_segment])
ax1.axvline(x=60,color='black',linestyle="--")
ax1.grid()
ax1.set_title('Condenser Air Outlet Pressure Distribution, N=60')
ax1.set_ylabel('P [kPa]')
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Loaction along tube (sgement) [-]')
fig.savefig('Performance_Pao.pdf')
#Q
fig, ax1 = plt.subplots()
ax1.plot(x,Qdist_dict[N_segment])
ax1.axvline(x=60,color='black',linestyle="--")
ax1.grid()
ax1.set_title('Condenser Capacity Distribution, N=60')
ax1.set_ylabel('Q [kW]')
# ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Loaction along tube (sgement) [-]')
fig.savefig('Performance_Q.pdf')
#htc_ref
#Q
fig, ax1 = plt.subplots()
ax1.plot(x,htcref_dist_dict[N_segment])
ax1.axvline(x=60,color='black',linestyle="--")
ax1.grid()
ax1.set_title('Condenser refrigerant HTC Distribution, N=60')
ax1.set_ylabel('heat transfer coefficient [kW/m^2-K]')
# ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax1.set_xlabel('Loaction along tube (sgement) [-]')
fig.savefig('Performance_htcref.pdf')
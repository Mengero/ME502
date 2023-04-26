import matplotlib.pyplot as plt
import numpy as np

N_elemptube = 40
T_ai = 27

# Data for Chen_1966 method
loaded_data = np.load('output_chen_1966.npz',allow_pickle=True)
T_ro_chen_1966 = loaded_data['T_ro']
P_ro_chen_1966 = loaded_data['P_ro']
x_ro_chen_1966 = loaded_data['x_ro']
T_ao_chen_1966 = loaded_data['T_ao']
P_ao_chen_1966 = loaded_data['P_ao']
Q_tot_chen_1966 = loaded_data['Q_tot']
T_ro_dict_chen_1966 = loaded_data['T_ro_dict'].item()
P_ro_dict_chen_1966 = loaded_data['P_ro_dict'].item()
x_ro_dict_chen_1966 = loaded_data['x_ro_dict'].item()
Q_dist_dict_chen_1966 = loaded_data['Q_dist_dict'].item()
T_r_dist_dict_chen_1966 = loaded_data['T_r_dist_dict'].item()
T_a_dist_dict_chen_1966 = loaded_data['T_a_dist_dict'].item()
P_r_dist_dict_chen_1966 = loaded_data['P_r_dist_dict'].item()
P_a_dist_dict_chen_1966 = loaded_data['P_a_dist_dict'].item()
DELTAP_air_dist_dict_chen_1966 = loaded_data['DELTAP_air_dist_dict'].item()
DELTAP_ref_dist_dict_chen_1966 = loaded_data['DELTAP_ref_dist_dict'].item()
x_r_dist_dict_chen_1966 = loaded_data['x_r_dist_dict'].item()
htc_ref_dist_dict_chen_1966 = loaded_data['htc_ref_dist_dict'].item()
Q_dist_chen_1966 = np.concatenate((np.flip(Q_dist_dict_chen_1966[1]),np.flip(Q_dist_dict_chen_1966[2]),np.flip(Q_dist_dict_chen_1966[3]),np.flip(Q_dist_dict_chen_1966[4])),axis=None)
T_r_dist_chen_1966 = np.concatenate((np.flip(T_r_dist_dict_chen_1966[1]),np.flip(T_r_dist_dict_chen_1966[2]),np.flip(T_r_dist_dict_chen_1966[3]),np.flip(T_r_dist_dict_chen_1966[4])),axis=None)
P_r_dist_chen_1966 = np.concatenate((np.flip(P_r_dist_dict_chen_1966[1]),np.flip(P_r_dist_dict_chen_1966[2]),np.flip(P_r_dist_dict_chen_1966[3]),np.flip(P_r_dist_dict_chen_1966[4])),axis=None)
T_ai_dist_chen_1966 = np.concatenate((T_a_dist_dict_chen_1966[2],np.flip(T_a_dist_dict_chen_1966[3]),T_a_dist_dict_chen_1966[4],np.flip(np.ones(N_elemptube)*T_ai)),axis=None)
T_ao_dist_chen_1966 = np.concatenate((T_a_dist_dict_chen_1966[1],np.flip(T_a_dist_dict_chen_1966[2]),T_a_dist_dict_chen_1966[3],np.flip(T_a_dist_dict_chen_1966[4])),axis=None)
x_r_dist_chen_1966 = np.concatenate((np.flip(x_r_dist_dict_chen_1966[1]),np.flip(x_r_dist_dict_chen_1966[2]),np.flip(x_r_dist_dict_chen_1966[3]),np.flip(x_r_dist_dict_chen_1966[4])),axis=None)
htc_ref_dist_chen_1966 = np.concatenate((np.flip(htc_ref_dist_dict_chen_1966[1]),np.flip(htc_ref_dist_dict_chen_1966[2]),np.flip(htc_ref_dist_dict_chen_1966[3]),np.flip(htc_ref_dist_dict_chen_1966[4])),axis=None)

# Data for Kim_2013 method
loaded_data = np.load('output_kim_2013.npz',allow_pickle=True)
T_ro_kim_2013 = loaded_data['T_ro']
P_ro_kim_2013 = loaded_data['P_ro']
x_ro_kim_2013 = loaded_data['x_ro']
T_ao_kim_2013 = loaded_data['T_ao']
P_ao_kim_2013 = loaded_data['P_ao']
Q_tot_kim_2013 = loaded_data['Q_tot']
T_ro_dict_kim_2013 = loaded_data['T_ro_dict'].item()
P_ro_dict_kim_2013 = loaded_data['P_ro_dict'].item()
x_ro_dict_kim_2013 = loaded_data['x_ro_dict'].item()
Q_dist_dict_kim_2013 = loaded_data['Q_dist_dict'].item()
T_r_dist_dict_kim_2013 = loaded_data['T_r_dist_dict'].item()
T_a_dist_dict_kim_2013 = loaded_data['T_a_dist_dict'].item()
P_r_dist_dict_kim_2013 = loaded_data['P_r_dist_dict'].item()
P_a_dist_dict_kim_2013 = loaded_data['P_a_dist_dict'].item()
DELTAP_air_dist_dict_kim_2013 = loaded_data['DELTAP_air_dist_dict'].item()
DELTAP_ref_dist_dict_kim_2013 = loaded_data['DELTAP_ref_dist_dict'].item()
x_r_dist_dict_kim_2013 = loaded_data['x_r_dist_dict'].item()
htc_ref_dist_dict_kim_2013 = loaded_data['htc_ref_dist_dict'].item()
Q_dist_kim_2013 = np.concatenate((np.flip(Q_dist_dict_kim_2013[1]),np.flip(Q_dist_dict_kim_2013[2]),np.flip(Q_dist_dict_kim_2013[3]),np.flip(Q_dist_dict_kim_2013[4])),axis=None)
T_r_dist_kim_2013 = np.concatenate((np.flip(T_r_dist_dict_kim_2013[1]),np.flip(T_r_dist_dict_kim_2013[2]),np.flip(T_r_dist_dict_kim_2013[3]),np.flip(T_r_dist_dict_kim_2013[4])),axis=None)
P_r_dist_kim_2013 = np.concatenate((np.flip(P_r_dist_dict_kim_2013[1]),np.flip(P_r_dist_dict_kim_2013[2]),np.flip(P_r_dist_dict_kim_2013[3]),np.flip(P_r_dist_dict_kim_2013[4])),axis=None)
T_ai_dist_kim_2013 = np.concatenate((T_a_dist_dict_kim_2013[2],np.flip(T_a_dist_dict_kim_2013[3]),T_a_dist_dict_kim_2013[4],np.flip(np.ones(N_elemptube)*T_ai)),axis=None)
T_ao_dist_kim_2013 = np.concatenate((T_a_dist_dict_kim_2013[1],np.flip(T_a_dist_dict_kim_2013[2]),T_a_dist_dict_kim_2013[3],np.flip(T_a_dist_dict_kim_2013[4])),axis=None)
x_r_dist_kim_2013 = np.concatenate((np.flip(x_r_dist_dict_kim_2013[1]),np.flip(x_r_dist_dict_kim_2013[2]),np.flip(x_r_dist_dict_kim_2013[3]),np.flip(x_r_dist_dict_kim_2013[4])),axis=None)
htc_ref_dist_kim_2013 = np.concatenate((np.flip(htc_ref_dist_dict_kim_2013[1]),np.flip(htc_ref_dist_dict_kim_2013[2]),np.flip(htc_ref_dist_dict_kim_2013[3]),np.flip(htc_ref_dist_dict_kim_2013[4])),axis=None)

x = np.array(range(1,4*N_elemptube+1))
# Plot T_r
fig, ax1 = plt.subplots()
ax1.plot(x,T_r_dist_kim_2013,label='Kim_2013')
ax1.plot(x,T_r_dist_chen_1966,label='Chen_1966')
ax1.axvline(x=40,color='black',linestyle="--")
ax1.axvline(x=80,color='black',linestyle="--")
ax1.axvline(x=120,color='black',linestyle="--")
ax1.set_xlim(0,160)
ax1.legend(loc='upper left')
ax1.grid()
ax1.set_title('Refrigerant Temperature along circuit, N_segment=40')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Nsegment_T_r.pdf')

# Plot T_ai
fig, ax1 = plt.subplots()
ax1.plot(x,T_ai_dist_kim_2013,label='Kim_2013')
ax1.plot(x,T_ai_dist_chen_1966,label='Chen_1966')
ax1.axvline(x=40,color='black',linestyle="--")
ax1.axvline(x=80,color='black',linestyle="--")
ax1.axvline(x=120,color='black',linestyle="--")
ax1.set_xlim(0,160)
ax1.legend(loc='upper left')
ax1.grid()
ax1.set_title('Air Inlet Temperature along circuit, N_segment=40')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Nsegment_T_ai.pdf')

# Plot T_ao
fig, ax1 = plt.subplots()
ax1.plot(x,T_ao_dist_kim_2013,label='Kim_2013')
ax1.plot(x,T_ao_dist_chen_1966,label='Chen_1966')
ax1.axvline(x=40,color='black',linestyle="--")
ax1.axvline(x=80,color='black',linestyle="--")
ax1.axvline(x=120,color='black',linestyle="--")
ax1.set_xlim(0,160)
ax1.legend(loc='upper left')
ax1.grid()
ax1.set_title('Air Outlet Temperature along circuit, N_segment=40')
ax1.set_ylabel('T [C]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Nsegment_T_ao.pdf')

# Plot P_r
fig, ax1 = plt.subplots()
ax1.plot(x,P_r_dist_kim_2013,label='Kim_2013')
ax1.plot(x,P_r_dist_chen_1966,label='Chen_1966')
ax1.axvline(x=40,color='black',linestyle="--")
ax1.axvline(x=80,color='black',linestyle="--")
ax1.axvline(x=120,color='black',linestyle="--")
ax1.set_xlim(0,160)
ax1.legend(loc='upper right')
ax1.grid()
ax1.set_title('Refrigerant Pressure along circuit, N_segment=40')
ax1.set_ylabel('P [kPa]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Nsegment_P_r.pdf')

# Plot x_r
fig, ax1 = plt.subplots()
ax1.plot(x,x_r_dist_kim_2013,label='Kim_2013')
ax1.plot(x,x_r_dist_chen_1966,label='Chen_1966')
ax1.axvline(x=40,color='black',linestyle="--")
ax1.axvline(x=80,color='black',linestyle="--")
ax1.axvline(x=120,color='black',linestyle="--")
ax1.set_xlim(0,160)
ax1.legend(loc='upper right')
ax1.grid()
ax1.set_title('Refrigerant Vapor Quality along circuit, N_segment=40')
ax1.set_ylabel('x [-]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Nsegment_x_r.pdf')


# # Plot HTC_r
fig, ax1 = plt.subplots()
ax1.plot(x,htc_ref_dist_kim_2013,label='Kim_2013')
ax1.plot(x,htc_ref_dist_chen_1966,label='Chen_1966')
ax1.axvline(x=40,color='black',linestyle="--")
ax1.axvline(x=80,color='black',linestyle="--")
ax1.axvline(x=120,color='black',linestyle="--")
ax1.set_xlim(0,160)
ax1.legend(loc='upper right')
ax1.grid()
ax1.set_title('Refrigerant heat transfer coefficient along circuit, N_segment=40')
ax1.set_ylabel('htc [kW/m^2-K]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Nsegment_htc_r.pdf')

# # Plot Q
fig, ax1 = plt.subplots()
ax1.plot(x,Q_dist_kim_2013,label='Kim_2013')
ax1.plot(x,Q_dist_chen_1966,label='Chen_1966')
ax1.axvline(x=40,color='black',linestyle="--")
ax1.axvline(x=80,color='black',linestyle="--")
ax1.axvline(x=120,color='black',linestyle="--")
ax1.set_xlim(0,160)
ax1.legend(loc='upper right')
ax1.grid()
ax1.set_title('Energy Transfer Rate along circuit, N_segment=40')
ax1.set_ylabel('Q [kW]')
ax1.set_xlabel('Circuit Location')
fig.savefig('Nsegment_Q.pdf')

import numpy as np
import CoolProp.CoolProp as CP
from simulation_utils.htc_2ph_r1234yf_Kim_2014 import htc_2ph_r1234yf_Kim_2014
from simulation_utils.htc_2ph_y1234yf_Shah import htc_2ph_y1234yf_shah
import matplotlib.pyplot as plt

def htc_vs(fld,x_r,G_ref,P_r,Dh_port,asp_ratio):

    '''
        Calculate htc based on different metods, variables: vapor quality
        mass flux G, and saturation temperature T_sat
    '''
    htc_1 = htc_2ph_r1234yf_Kim_2014(fld,G_ref,Dh_port,asp_ratio,P_r,x_r)       # [W/m^2-K]
    htc_2 = htc_2ph_y1234yf_shah(fld,G_ref,Dh_port,P_r,x_r)                # [W/m^2-K]

    return htc_1, htc_2


T_ri = 40
P_r = 10.17*101.325     # [kPa]
G_ref = 200


Dh_port = 0.96e-3
asp_ratio = 1
fld = 'R1234yf'
T_sat = CP.PropsSI('T','P',P_r*1e3,'Q',0.5,fld)-273.15
print(T_sat)
# calculate for quality
htc_kim_Q=np.array([]); htc_shah_Q=np.array([])
x_r = np.linspace(0.1,0.99,num=30)
for x_r_tmp in x_r:
    [htc_kim_Q_tmp,htc_shah_Q_tmp] = htc_vs(fld,x_r_tmp,G_ref,P_r,Dh_port,asp_ratio)
    htc_kim_Q = np.append(htc_kim_Q,htc_kim_Q_tmp); htc_shah_Q = np.append(htc_shah_Q,htc_shah_Q_tmp)
fig, ax = plt.subplots()
ax.plot(x_r,htc_kim_Q,label="Kim (2014)")
ax.plot(x_r,htc_shah_Q,label='Shah (2019)')
ax.set_xlabel('Vapor Quality [-]')
ax.set_ylabel('Heat Transfer Coefficient [W/m^2-K]')
ax.set_title('T_sat=40.46[C], G=200[kg/m^2-s], R1234yf')
legend = ax.legend(loc='upper left')
legend.get_frame().set_facecolor('C0')
fig.savefig('htc_Q_comp.pdf')


# T_sat
x_r = 0.5
P_r = np.linspace(1200,1600,num=100)       # [kPa]
htc_kim_Tsat=np.array([]); htc_shah_Tsat=np.array([])
for P_r_tmp in P_r:
    [htc_kim_Tsat_tmp,htc_shah_Tsat_tmp] = htc_vs(fld,x_r,G_ref,P_r_tmp,Dh_port,asp_ratio)
    htc_kim_Tsat = np.append(htc_kim_Tsat,htc_kim_Tsat_tmp); htc_shah_Tsat = np.append(htc_shah_Tsat,htc_shah_Tsat_tmp)
T_sat = CP.PropsSI('T','P',P_r*1e3,'Q',0.5,fld)-273.15
fig, ax = plt.subplots()
ax.plot(T_sat,htc_kim_Tsat,label="Kim (2014)")
ax.plot(T_sat,htc_shah_Tsat,label='Shah (2019)')
ax.set_xlabel('Saturation Temperature [C]')
ax.set_ylabel('Heat Transfer Coefficient [W/m^2-K]')
ax.set_title('x=0.5, G=200[kg/m^2-s], R1234yf')
legend = ax.legend(loc='upper right')
legend.get_frame().set_facecolor('C0')
fig.savefig('htc_Tsat_comp.pdf')

# G
x_r = 0.5
G_ref = np.linspace(53,1403,num=100) # [kg/m^2-s]
P_r = 10.17*101.325       # [kPa]
htc_kim_G=np.array([]); htc_shah_G=np.array([])
for G_ref_tmp in G_ref:
    [htc_kim_Gref_tmp,htc_shah_Gref_tmp] = htc_vs(fld,x_r,G_ref_tmp,P_r,Dh_port,asp_ratio)
    htc_kim_G = np.append(htc_kim_G,htc_kim_Gref_tmp); htc_shah_G = np.append(htc_shah_G,htc_shah_Gref_tmp)
T_sat = CP.PropsSI('T','P',P_r*1e3,'Q',0.5,fld)-273.15
fig, ax = plt.subplots()
ax.plot(G_ref,htc_kim_G,label="Kim (2014)")
ax.plot(G_ref,htc_shah_G,label='Shah (2019)')
ax.set_xlabel('Mass Flux [kg/m^2-s]')
ax.set_ylabel('Heat Transfer Coefficient [W/m^2-K]')
ax.set_title('T_r=40.46[C], x=0.5, R1234yf')
legend = ax.legend(loc='upper left')
legend.get_frame().set_facecolor('C0')
fig.savefig('htc_G_comp.pdf')


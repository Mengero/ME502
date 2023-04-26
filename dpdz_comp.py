import numpy as np
import CoolProp.CoolProp as CP
from simulation_utils.dpdz_f_2ph_Del_Col_2013 import dpdz_f_2ph_Del_Col_2013
from simulation_utils.dpdz_f_2ph_Lopez_2014 import dpdz_f_2ph_Lopez_2014
import matplotlib.pyplot as plt

def DP_vs(fld,roughness,Dh,G_ref,P_r,x_r,L):

    rho_ref = CP.PropsSI('D','P',P_r*1e3,'Q',x_r,fld)
    T_sat=CP.PropsSI('T','P',P_r*1e3,'Q',x_r,fld)-273.15
    dpdz_Del = dpdz_f_2ph_Del_Col_2013(fld,roughness,Dh,G_ref,P_r,x_r)
    dpdz_Lopez = dpdz_f_2ph_Lopez_2014(fld,G_ref,T_sat,x_r,Dh)

    DP_static_Del = dpdz_Del/Dh * G_ref**2/(2*rho_ref)
    DP_static_Lopez = dpdz_Lopez

    return DP_static_Del,DP_static_Lopez

G_ref = 200         # [kg/m^2-s]
fld = 'R1234yf'
T_sat = 36          # [C]
Dh = 1.23e-3       # [m]
L = 0.249       # [m]
roughness=1e-6

# Quality Comparison
P_r = CP.PropsSI('P','T',T_sat+273.15,'Q',0.5,fld)/1e3      # [kPa]
DP_Del_Q=np.array([]); DP_Lopez_Q=np.array([])
x_r = np.linspace(0.1,0.99,num=30)
for x_r_tmp in x_r:
    [DP_Del_Q_tmp,DP_Lopez_Q_tmp] = DP_vs(fld,roughness,Dh,G_ref,P_r,x_r_tmp,L)
    DP_Del_Q = np.append(DP_Del_Q,DP_Del_Q_tmp/1e3); DP_Lopez_Q = np.append(DP_Lopez_Q,DP_Lopez_Q_tmp/1e3)
fig, ax = plt.subplots()
ax.plot(x_r,DP_Del_Q,label="Del")
ax.plot(x_r,DP_Lopez_Q,label='Lopez')
ax.set_xlabel('Vapor Quality [-]')
ax.set_ylabel('Pressure Drop [kPa/m]')
ax.set_title('G=200[kg/m^2-s], T_sat=36[C], R1234yf')
legend = ax.legend(loc='upper left')
legend.get_frame().set_facecolor('C0')
fig.savefig('DP_Q_comp.pdf')

# G_ref comparison
P_r = CP.PropsSI('P','T',T_sat+273.15,'Q',0.5,fld)/1e3      # [kPa]
DP_Del_Q=np.array([]); DP_Lopez_Q=np.array([])
x_r = 0.5
G_ref = np.linspace(100,1400,num=30)
for G_ref_tmp in G_ref:
    [DP_Del_Q_tmp,DP_Lopez_Q_tmp] = DP_vs(fld,roughness,Dh,G_ref_tmp,P_r,x_r,L)
    DP_Del_Q = np.append(DP_Del_Q,DP_Del_Q_tmp/1e3); DP_Lopez_Q = np.append(DP_Lopez_Q,DP_Lopez_Q_tmp/1e3)
fig, ax = plt.subplots()
ax.plot(G_ref,DP_Del_Q,label="Del")
ax.plot(G_ref,DP_Lopez_Q,label='Lopez')
ax.set_xlabel('Mass Velocity [mg/m^2-s]')
ax.set_ylabel('Pressure Drop [kPa/m]')
ax.set_title('x=0.5, T_sat=36[C], R1234yf')
legend = ax.legend(loc='upper left')
legend.get_frame().set_facecolor('C0')
fig.savefig('DP_G_comp.pdf')

# T_sat
T_sat = np.linspace(20,50,num=30)
P_r = CP.PropsSI('P','T',T_sat+273.15,'Q',0.5,fld)/1e3      # [kPa]
DP_Del_Q=np.array([]); DP_Lopez_Q=np.array([])
x_r = 0.5
G_ref = 200
for P_r_tmp in P_r:
    [DP_Del_Q_tmp,DP_Lopez_Q_tmp] = DP_vs(fld,roughness,Dh,G_ref,P_r_tmp,x_r,L)
    DP_Del_Q = np.append(DP_Del_Q,DP_Del_Q_tmp/1e3); DP_Lopez_Q = np.append(DP_Lopez_Q,DP_Lopez_Q_tmp/1e3)
fig, ax = plt.subplots()
ax.plot(T_sat,DP_Del_Q,label="Del")
ax.plot(T_sat,DP_Lopez_Q,label='Lopez')
ax.set_xlabel('Saturation Temperature [C]')
ax.set_ylabel('Pressure Drop [kPa/m]')
ax.set_title('G=200[kg/m^2-s], x=0.5, R1234yf')
legend = ax.legend(loc='upper right')
legend.get_frame().set_facecolor('C0')
fig.savefig('DP_Tsat_comp.pdf')
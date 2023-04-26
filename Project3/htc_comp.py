from simulation_utils.htc_2ph_Kim_2013 import htc_2ph_Kim_2013
from simulation_utils.htc_2ph_chen_1966 import htc_2ph_chen_1966
import numpy as np
from evaporator import evaporator
import CoolProp.CoolProp as CP
from simulation_utils.airsidehtc_chang_wang import airsidehtc_chang_wang
import matplotlib.pyplot as plt

Ref = 'R1234yf'
x_ri = 0.15
T_sat_ri = 7
m_dot_ref = 25e-3
T_ai = 27
P_ai = 99.5
vol_dot_air = 500
N_segment = 40
G=200
self = evaporator(x_ri,T_sat_ri,m_dot_ref,T_ai,P_ai,vol_dot_air,N_segment)
[htc_air, eta_f, eta_t] = airsidehtc_chang_wang(P_ai,self.k_w,self.theta_louver,self.P_fin,self.P_louver,\
    self.h_fin,self.D_fin,self.L_louver,self.D_tube,self.t_tube,self.t_fin,self.T_ai,self.vel_air)
htc_air = htc_air/1e3
R_air = 1 / (self.htc_air * self.A_air_elem * eta_t)

def htc_vs(self,Ref,x_r,G,P_r,Dh,AspRat):
    T_r = CP.PropsSI('T','P',P_r*1e3,'Q',x_r,Ref)-273.15
    h_ri = CP.PropsSI('H','P',P_r*1e3,'Q',x_r,Ref)
    htc_1,_ = htc_2ph_chen_1966(Ref,T_r,P_r,x_r,G,Dh,self.R_air,T_ai,self.A_ref_elem,AspRat)              # [kW/m^2-K]
    htc_1 = htc_1/1e3
    q_H = 2.4e4      # [W/m^2]
    htc_2 = htc_2ph_Kim_2013(Ref,q_H,G,T_r,x_r,Dh,self.P_tube_outer,self.P_port_tot,AspRat)/1e3

    return htc_1, htc_2

# --------------------- x_r ---------------------- #
htc_kim_Q=np.array([]); htc_chen_Q=np.array([])
x_r = np.linspace(0.1,0.99,num=30)
for x_r_tmp in x_r:
    [htc_chen_Q_tmp,htc_kim_Q_tmp] = htc_vs(self,Ref,x_r_tmp,self.G_ref_elem,self.P_ri,self.Dh_port,self.asp_ratio)
    htc_kim_Q = np.append(htc_kim_Q,htc_kim_Q_tmp); htc_chen_Q = np.append(htc_chen_Q,htc_chen_Q_tmp)
fig, ax = plt.subplots()
ax.plot(x_r,htc_kim_Q,label="Kim (2013)")
ax.plot(x_r,htc_chen_Q,label='Chen (1966)')
ax.set_xlabel('Vapor Quality [-]')
ax.set_ylabel('Heat Transfer Coefficient [kW/m^2-K]')
ax.set_title('T_sat=7[C], G=143[kg/m^2-s], R1234yf')
legend = ax.legend(loc='upper right')
legend.get_frame().set_facecolor('C0')
fig.savefig('htc_Q_comp.pdf')

# --------------------- G ----------------------- #
G_ref = np.linspace(53,1403,num=100) # [kg/m^2-s]
htc_kim_G=np.array([]); htc_chen_G=np.array([])
for G_ref_tmp in G_ref:
    [htc_chen_G_tmp,htc_kim_G_tmp] = htc_vs(self,Ref,x_ri,G_ref_tmp,self.P_ri,self.Dh_port,self.asp_ratio)
    htc_kim_G = np.append(htc_kim_G,htc_kim_G_tmp); htc_chen_G = np.append(htc_chen_G,htc_chen_G_tmp)
fig, ax = plt.subplots()
ax.plot(G_ref,htc_kim_G,label="Kim (2013)")
ax.plot(G_ref,htc_chen_G,label='Chen (1966)')
ax.set_xlabel('Mass Flux [kg/m^2-s]')
ax.set_ylabel('Heat Transfer Coefficient [kW/m^2-K]')
ax.set_title('T_sat=7[C], x=0.15, R1234yf')
legend = ax.legend(loc='upper left')
legend.get_frame().set_facecolor('C0')
fig.savefig('htc_G_comp.pdf')

# --------------------- T_sat ----------------------- #
P_r = np.linspace(200,500,num=100)       # [kPa]
htc_kim_Tsat=np.array([]); htc_chen_Tsat=np.array([])
for P_r_tmp in P_r:
    [htc_chen_Tsat_tmp,htc_kim_Tsat_tmp] = htc_vs(self,Ref,x_ri,self.G_ref_elem,P_r_tmp,self.Dh_port,self.asp_ratio)
    htc_kim_Tsat = np.append(htc_kim_Tsat,htc_kim_Tsat_tmp); htc_chen_Tsat = np.append(htc_chen_Tsat,htc_chen_Tsat_tmp)
fig, ax = plt.subplots()
T_sat = CP.PropsSI('T','P',P_r*1e3,'Q',x_ri,Ref)-273.15
ax.plot(P_r,htc_kim_Tsat,label="Kim (2013)")
ax.plot(P_r,htc_chen_Tsat,label='Chen (1966)')
ax.set_xlabel('Saturation Temperature [C]')
ax.set_ylabel('Heat Transfer Coefficient [kW/m^2-K]')
ax.set_title('x=0.15, G=143[kg/m^2-s], R1234yf')
legend = ax.legend(loc='upper left')
legend.get_frame().set_facecolor('C0')
fig.savefig('htc_Tsat_comp.pdf')
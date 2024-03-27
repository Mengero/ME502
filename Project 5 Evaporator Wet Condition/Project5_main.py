from evaporator import evaporator_system
import pickle
import numpy as np
'''
    The impact of air inlet relative humidity and volumetric flow rate is investigated.
        RH_ai: 0.2 ~ 0.8
'''
# --------------------- Base Line --------------------- #
T_cro = 40          # [C]
P_cro = 1250        # [kPa]
m_dot_ref = 35      # [g/s]
Tsh_ro = 8          # [C]
T_ai = 35           # [C]
P_ai = 99.5         # [kPa]
RH_ai = 0.4
m_dot_air = 9        # [kg/min] 
N_segment = 40
Ref = 'R1234yf'

N_slab = 10
result_baseline = evaporator_system(Ref,T_cro,P_cro,m_dot_ref,Tsh_ro,T_ai,P_ai,RH_ai,m_dot_air,N_segment)

with open("baseline_performance.pkl", "wb") as f:
    pickle.dump((result_baseline), f)

# -------------------- RH_ai Impact -------------------- #
result_RH_ai = {}
T_cro = 40          # [C]
P_cro = 1250        # [kPa]
m_dot_ref = 35      # [g/s]
Tsh_ro = 8          # [C]
T_ai = 35           # [C]
P_ai = 99.5         # [kPa]
RH_ai = 0.4
m_dot_air = 9        # [kg/min] 
Ref = 'R1234yf'
N_segment = 40
for RH_ai in np.linspace(0.2,0.8,num=10):
    print('solving cases for RH_ai = {}'.format(RH_ai))
    result_tmp = evaporator_system(Ref,T_cro,P_cro,m_dot_ref,Tsh_ro,T_ai,P_ai,RH_ai,m_dot_air,N_segment)
    result_RH_ai[RH_ai] = result_tmp

with open("RH_ai_impact.pkl", "wb") as f:
    pickle.dump((result_RH_ai), f)

print('Simulation Finished')
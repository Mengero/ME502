from cooling_coil import cooling_coil
from cooling_coil import cc_system
from simulation_utils.airsidehtc_chang_wang import airsidehtc_chang_wang
import pickle
import numpy as np
'''
    The impact of air inlet relative humidity and volumetric flow rate is investigated.
        RH_ai: 0.2 ~ 0.8
        vol_dot_a: 250 ~ 600 [CFM]
'''
# -------------------- Overall Performance ---------------------- #
T_ri = 0                        # [C]
P_ri = 150                      # [kPa]
vol_dot_r = 30                  # [L/min]
T_ai = 27                       # [C]
P_ai = 99.5                     # [kPa]
RH_ai = 0.5                    # [-], relative humidity
vol_dot_a = 500               # [CFM]
N_segment = 40
self=cooling_coil(T_ri,P_ri,vol_dot_r,T_ai,P_ai,vol_dot_a,RH_ai,N_segment)
result_overall = cc_system(self,T_ri,P_ri,T_ai,P_ai,RH_ai)

with open("overall_performance.pkl", "wb") as f:
    pickle.dump((result_overall), f)

# -------------------- RH_ai impact evaluation ------------------- #
result_RH_ai = {}

T_ri = 0                        # [C]
P_ri = 150                      # [kPa]
vol_dot_r = 30                  # [L/min]
T_ai = 27                       # [C]
P_ai = 99.5                     # [kPa]
RH_ai = 0.5                    # [-], relative humidity
vol_dot_a = 500               # [CFM]
N_segment = 40
for RH_ai in np.linspace(0.2,0.8,num=10):
    print('solving cases for RH_ai = {}'.format(RH_ai))
    self=cooling_coil(T_ri,P_ri,vol_dot_r,T_ai,P_ai,vol_dot_a,RH_ai,N_segment)
    result_tmp = cc_system(self,T_ri,P_ri,T_ai,P_ai,RH_ai)
    result_RH_ai[RH_ai] = result_tmp

with open("RH_ai_impact.pkl", "wb") as f:
    pickle.dump((result_RH_ai), f)

# -------------------- vol_dot_a impact evaluation ------------------- #
result_vol_dot_a = {}

T_ri = 0                        # [C]
P_ri = 150                      # [kPa]
vol_dot_r = 30                  # [L/min]
T_ai = 27                       # [C]
P_ai = 99.5                     # [kPa]
RH_ai = 0.5                    # [-], relative humidity
vol_dot_a = 500               # [CFM]
N_segment = 40
for vol_dot_a in np.linspace(250,600,num=10):
    print('solving cases for vol_dot_a = {}'.format(vol_dot_a))
    self=cooling_coil(T_ri,P_ri,vol_dot_r,T_ai,P_ai,vol_dot_a,RH_ai,N_segment)
    result_tmp = cc_system(self,T_ri,P_ri,T_ai,P_ai,RH_ai)
    result_vol_dot_a[vol_dot_a] = result_tmp

# -------------------- Save variables to file ---------------------- #
with open("vol_dot_a_impact.pkl", "wb") as f:
    pickle.dump((result_vol_dot_a), f)
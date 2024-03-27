from BESS_geo import BESS
from BESS_control import BESS_control
import numpy as np


def BESS_main(Q_external,Q_batterygenerating):
    BESS_geo = BESS(Q_external,Q_batterygenerating)
    t = np.array([]); Tdb_air = np.array([]); Tdp_air = np.array([])
    T_eg = np.array([]); T_bat = np.array([]); T_dev = np.array([])
    T_frame = np.array([]); T_wall = np.array([])
    Q_battery = np.array([]); Q_air = np.array([]); Q_wall = np.array([])
    Q_devices = np.array([]); Q_frame = np.array([]); Q_EG = np.array([])
    Q_ht_air_frame = np.array([]); Q_ht_eg_battery = np.array([])
    Q_ht_air_wall = np.array([]); Q_ht_air_device = np.array([]); Q_ht_air_battery = np.array([])

    Time = 0                # [hr]
    DTime = 2               # [s]


    Tdb_air_ini = BESS_geo.Tdb_air_0; Tdp_air_ini = BESS_geo.Tdp_air_0; T_eg_ini = BESS_geo.T_eg_0
    T_bat_ini = BESS_geo.T_bat_0; T_dev_ini = BESS_geo.T_dev_0; T_frame_ini = BESS_geo.T_frame_0
    T_wall_ini = BESS_geo.T_wall_0
    air_cooling = 0
    bat_cooling = 0
    while Time <= 20:
        Time += DTime/3600
        t = np.append(t,Time)
        if Tdb_air_ini>=BESS_geo.Tdb_air_max or Tdp_air_ini>=BESS_geo.Tdp_air_max:
            air_cooling = 1
        elif Tdb_air_ini<=BESS_geo.Tdb_air_min and Tdp_air_ini<=BESS_geo.Tdp_air_min:
            air_cooling = 0

        if T_bat_ini>=BESS_geo.T_bat_max:
            bat_cooling = 1
        elif T_bat_ini<=BESS_geo.T_bat_max-0.5:
            bat_cooling = 0

        [T_bat_ini,T_eg_ini,Tdb_air_ini,Tdp_air_ini,T_wall_ini,T_dev_ini,T_frame_ini,\
         Q_battery_tmp,Q_air_tmp,Q_wall_tmp,Q_devices_tmp,Q_frame_tmp,Q_EG_tmp,\
         Q_ht_air_battery_tmp,Q_ht_eg_battery_tmp,Q_ht_air_wall_tmp,Q_ht_air_device_tmp,Q_ht_air_frame_tmp] = BESS_control(\
         BESS_geo,Tdb_air_ini,Tdp_air_ini,T_eg_ini,T_bat_ini,T_dev_ini,T_frame_ini,T_wall_ini,DTime,air_cooling,bat_cooling)
        Tdb_air=np.append(Tdb_air,Tdb_air_ini); Tdp_air=np.append(Tdp_air,Tdp_air_ini); T_eg=np.append(T_eg,T_eg_ini)
        T_bat=np.append(T_bat,T_bat_ini); T_dev=np.append(T_dev,T_dev_ini); T_frame=np.append(T_frame,T_frame_ini)
        T_wall=np.append(T_wall,T_wall_ini); Q_battery=np.append(Q_battery,Q_battery_tmp); Q_air=np.append(Q_air,Q_air_tmp)
        Q_wall=np.append(Q_wall,Q_wall_tmp); Q_devices=np.append(Q_devices,Q_devices_tmp); Q_frame=np.append(Q_frame,Q_frame_tmp)
        Q_EG=np.append(Q_EG,Q_EG_tmp); Q_ht_air_battery=np.append(Q_ht_air_battery,Q_ht_air_battery_tmp); Q_ht_eg_battery=np.append(Q_ht_eg_battery,Q_ht_eg_battery_tmp)
        Q_ht_air_wall=np.append(Q_ht_air_wall,Q_ht_air_wall_tmp);Q_ht_air_device=np.append(Q_ht_air_device,Q_ht_air_device_tmp)
        Q_ht_air_frame=np.append(Q_ht_air_frame,Q_ht_air_frame_tmp)

    print('finished')
    return t,Tdb_air,Tdp_air,T_eg,T_bat,T_dev,T_frame,T_wall,Q_battery,Q_air,Q_wall,Q_devices,Q_frame,Q_EG,\
           Q_ht_air_battery,Q_ht_eg_battery,Q_ht_air_wall,Q_ht_air_device,Q_ht_air_frame

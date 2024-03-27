import CoolProp.CoolProp as CP

def BESS_control(self,Tdb_air_ini,Tdp_air_ini,T_eg_ini,T_bat_ini,T_dev_ini,T_frame_ini,T_wall_ini,DTime,air_cooling,bat_cooling):
    '''
        Assume air pressure is constant in this control volume
    '''
    SHR = -0.019*Tdp_air_ini + 1.12             # [C], sensible heat ratio
    SHR = max(0,min(1,SHR))

    Q_batterygenerating = self.Q_batterygenerating
    Q_ht_air_battery = (Tdb_air_ini - T_bat_ini)/self.R_air_bat
    Q_ht_eg_battery = (T_eg_ini - T_bat_ini)/self.R_eg_bat
    Q_ht_air_wall = (Tdb_air_ini - T_wall_ini)/self.R_air_wall
    Q_ht_air_device = (Tdb_air_ini - T_dev_ini)/self.R_air_dev
    Q_ht_air_frame = (Tdb_air_ini - T_frame_ini)/self.R_air_frame
    Q_external = self.Q_external
    Q_devicegenerating = self.Power_device

    # Determine Control Strategy
    if bat_cooling==1:
        if air_cooling == 1:
            if (Tdb_air_ini>=self.Tdb_air_max) or (Tdp_air_ini>=self.Tdp_air_max):
                # Battery and air needs cooling, both evaporator is working
                Q_air_cooling = self.Q_air_cooling; Q_air_heater = 0
                print('Battery Cooling, Air Cooling and Dehumidifying, 2 Evaporator ON, # MODE 3 #')
            elif Tdp_air_ini<=self.Tdp_air_min and Tdb_air_ini>=self.Tdb_air_min:
                Q_air_cooling = self.Q_air_cooling; Q_air_heater = 0
                print('Battery Cooling, Air Cooling and Dehumidifying, 2 Evaporator ON, # MODE 3 #')
            elif Tdp_air_ini>=self.Tdp_air_min and Tdb_air_ini<=self.Tdb_air_min:
                # Only battery needs cooling, air needs heating, only first evaporator is working
                Q_air_cooling = self.Q_air_cooling; Q_air_heater = self.Q_air_heater
                print('Battery Cooling, Air Heating and Dehumidifying, 2 Evaporator ON, # MODE 3 #')
            else:
                # Only battery needs cooling, air needs cooling, only first evaporator is working
                Q_air_cooling = self.Q_air_cooling; Q_air_heater = 0
                print('Battery Cooling, Air Heating and Dehumidifying, 2 Evaporator ON, # MODE 3 #')
        else:
            Q_air_cooling = 0; Q_air_heater = 0
            print('Battery Cooling, Air Chilling, 1 Evaporator ON, # MODE 2 #')
        Q_eg_cooling = self.Q_cooling - Q_air_cooling; Q_eg_heater = 0
    else:
        if air_cooling == 1:
            if Tdb_air_ini>=self.Tdb_air_max or Tdp_air_ini>=self.Tdp_air_max:
                # Battery doesn't need cooling, but air needs cooling, both evaporator working but heater is on for EG
                Q_air_cooling = self.Q_air_cooling; Q_air_heater = 0
                Q_eg_cooling = self.Q_cooling - Q_air_cooling
                print('Battery Chilling, Air Cooling and Dehumidifying, 2 Evaporator ON, # MODE 3 #')
            elif Tdb_air_ini<=self.Tdb_air_min and Tdp_air_ini>=self.Tdp_air_min:
                # Battery doesn't needs cooling, but air needs heating
                Q_air_cooling = self.Q_air_cooling; Q_air_heater = self.Q_air_heater
                Q_eg_cooling = self.Q_cooling - Q_air_cooling
                print('Battery Chilling, Air heating and Dehumidifying, 2 Evaporator ON, # MODE 3 #')
            elif Tdb_air_ini>=self.Tdb_air_min and Tdp_air_ini<=self.Tdp_air_min:
                # Battery doesn't needs cooling, but air needs heating
                Q_air_cooling = self.Q_air_cooling; Q_air_heater = 0
                Q_eg_cooling = self.Q_cooling - Q_air_cooling
                print('Battery Chilling, Air Cooling, 2 Evaporator ON, # MODE 3 #')
            else:
                Q_air_cooling = self.Q_air_cooling; Q_air_heater = 0
                Q_eg_cooling = self.Q_cooling - Q_air_cooling
                print('Battery Chilling, Air Cooling and Dehumidifying, 2 Evaporator ON, # MODE 3 #')

        else:
            Q_air_cooling = 0; Q_air_heater = 0
            Q_eg_cooling = 0
            print('Battery Chilling, Air Chilling, 0 Evaporator ON, # MODE 1 #')
    
    if T_bat_ini<=self.T_bat_min:
        Q_eg_heater = self.Q_eg_heater
    else:
        Q_eg_heater = 0


    # ----------------- Battery ------------------ #
    Q_battery = Q_batterygenerating + Q_ht_air_battery + Q_ht_eg_battery
    T_bat_out = T_bat_ini + Q_battery*DTime/(self.cp_bat*self.m_bat)
    # ----------------- Water_EG ------------------ #
    Q_EG = Q_eg_heater-Q_eg_cooling-Q_ht_eg_battery + 0.8*self.Power_pump
    T_eg_out = T_eg_ini + Q_EG*DTime/(self.cp_eg*self.m_eg)
    # ----------------- air ------------------ #
    Q_air_lat = Q_air_cooling*(1-SHR)
    Q_air_sen = Q_air_cooling*SHR
    Q_air = - Q_ht_air_wall - Q_ht_air_device - Q_ht_air_frame - Q_air_sen + Q_air_heater - Q_ht_air_battery + self.Power_fan + 0.2*self.Power_pump
    h_f = CP.PropsSI('H','T',Tdp_air_ini+273.15,'Q',0,'Water')/1e3; h_g = CP.PropsSI('H','T',Tdp_air_ini+273.15,'Q',1,'Water')/1e3
    h_fg_water = h_g - h_f
    m_dot_water = -Q_air_lat/h_fg_water*DTime + self.m_dot_water_in*DTime/1e3
    Vha = CP.HAPropsSI('Vha','Tdb',Tdb_air_ini+273.15,'Tdp',Tdp_air_ini+273.15,'P',self.P_air_0*1e3)
    Vda = CP.HAPropsSI('Vda','Tdb',Tdb_air_ini+273.15,'Tdp',Tdp_air_ini+273.15,'P',self.P_air_0*1e3)
    omega_ini = CP.HAPropsSI('Omega','Tdb',Tdb_air_ini+273.15,'Tdp',Tdp_air_ini+273.15,'P',self.P_air_0*1e3)
    m_air = self.V_air/Vha
    rho_H2O_ini = omega_ini/Vda
    m_water_ini = self.V_air*rho_H2O_ini
    m_da = m_air-m_water_ini
    omega_out = (m_water_ini + m_dot_water)/m_da
    cp_air = CP.HAPropsSI('C','Tdb',Tdb_air_ini+273.15,'Tdp',Tdp_air_ini+273.15,'P',self.P_air_0*1e3)/1e3
    Tdb_air_out = Tdb_air_ini + Q_air*DTime/(cp_air*m_da)
    Tdp_air_out = CP.HAPropsSI('Tdp','Omega',omega_out,'Tdb',Tdb_air_out+273.15,'P',self.P_air_0*1e3)-273.15
    # ----------------- Wall ------------------ #
    Q_wall = Q_external + Q_ht_air_wall
    T_wall_out = T_wall_ini + Q_wall*DTime/(self.cp_wall*self.m_wall)
    # ----------------- Devices ------------------ #
    Q_devices = Q_devicegenerating + Q_ht_air_device
    T_dev_out = T_dev_ini + Q_devices*DTime/(self.cp_dev*self.m_dev)
    # ----------------- Frame ------------------ #
    Q_frame = Q_ht_air_frame
    T_frame_out = T_frame_ini + Q_frame*DTime/(self.cp_frame*self.m_frame)

    return T_bat_out,T_eg_out,Tdb_air_out,Tdp_air_out,T_wall_out,T_dev_out,T_frame_out,Q_battery,Q_air,Q_wall,Q_devices,Q_frame,Q_EG,\
           Q_ht_air_battery,Q_ht_eg_battery,Q_ht_air_wall,Q_ht_air_device,Q_ht_air_frame
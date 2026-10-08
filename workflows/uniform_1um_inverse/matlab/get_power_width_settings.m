function settings = get_power_width_settings()
%GET_POWER_WIDTH_SETTINGS Shared optimization and width conversion settings.

settings.Initial_Voltage_V = 30;
settings.Max_Contrast = 4;
settings.Rho_Elec_Ohm_m = 1 / (75272 / 6);
settings.Device_Thickness_m = 220e-9;
settings.Ambient_Temp_K = 293.15;

end

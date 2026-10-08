function settings = get_power_width_settings()
%GET_POWER_WIDTH_SETTINGS Shared optimization and width conversion settings.

settings.Initial_Voltage_V = 30;
settings.Max_Contrast = 4;
settings.Rho_Elec_Ohm_m = 1 / (75272 / 6);
settings.Device_Thickness_m = 220e-9;
settings.Ambient_Temp_K = 293.15;

% In X mode, penalize second differences in interior-column and top/bottom
% edge-group power densities. Linear ramps have zero second difference;
% alternating-column oscillations are penalized.
% Columns 1 and Nx are excluded to retain independent edge compensation.
% The augmented least-squares system [C; lambda*D] x ~ [d; 0] minimizes
% ||C x - d||^2 + lambda^2 ||D x||^2. Set lambda to zero to disable this term.
settings.X_Smoothing_Lambda = 1e4;

end

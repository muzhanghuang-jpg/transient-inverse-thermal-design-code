function settings = get_thermal_target_settings()
%GET_THERMAL_TARGET_SETTINGS Central target definition for optimization steps.

    settings.Base_Temp_K = 900;
    settings.Y_Gradient_K_per_um = 2;
    % GradientDirection selects the ramp axis ('X' or 'Y') for the target,
    % group masks, and response-matrix construction.
    % The Y_Gradient_K_per_um field name is kept for backward compatibility
    % and represents the gradient magnitude along the chosen direction.
    settings.GradientDirection = 'X';
    settings.Target_Time_s = 13E-6;
    settings.Time_Step_s = 1E-6;
end

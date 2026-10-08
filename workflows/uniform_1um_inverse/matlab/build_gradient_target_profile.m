function [Target_Temp_ROI, Target_Temp_full, Y_um_full, y_ref_um] = ...
    build_gradient_target_profile(W_arr, H_arr, monitor_indices, baseTempK, gradientKPerUm)
%BUILD_GRADIENT_TARGET_PROFILE Build a Y-gradient target centered on the ROI.
%
% Positive gradient means larger Y coordinates have higher target temperature.

    Nx = numel(W_arr);
    Ny = numel(H_arr);

    Y_bl = zeros(1, Ny);
    Y_bl(1) = -H_arr(1);
    for j = 2:Ny
        Y_bl(j) = Y_bl(j-1) + H_arr(j-1);
    end

    yCenterUm = (Y_bl + H_arr / 2) / 1000;
    Y_um_full = repmat(yCenterUm, Nx, 1);

    roiYUm = Y_um_full(monitor_indices);
    y_ref_um = mean(roiYUm);

    Target_Temp_full = baseTempK + gradientKPerUm * (Y_um_full - y_ref_um);
    Target_Temp_ROI = Target_Temp_full(monitor_indices);
end

function [Target_Temp_ROI, Target_Temp_full, Coord_um_full, ref_um] = ...
    build_gradient_target_profile(W_arr, H_arr, monitor_indices, ...
    baseTempK, gradientKPerUm, gradientDirection)
%BUILD_GRADIENT_TARGET_PROFILE Build a 1D temperature ramp target.
%
% gradientDirection: 'X' or 'Y' (defaults to 'Y' for backward compatibility).
%   'X' makes the ramp vary along X (column index varies, target is constant
%   within each X-column); 'Y' makes the ramp vary along Y (row index varies,
%   target is constant within each Y-row).
%
% Coord_um_full is the coordinate along the selected ramp axis, X or Y.
% Some callers store this output as Y_um_full; the variable name does not
% determine the physical direction.

    if nargin < 6 || isempty(gradientDirection)
        gradientDirection = 'Y';
    end

    Nx = numel(W_arr);
    Ny = numel(H_arr);

    if strcmpi(gradientDirection, 'X')
        X_bl = zeros(1, Nx);
        X_bl(1) = -W_arr(1);
        for i = 2:Nx
            X_bl(i) = X_bl(i-1) + W_arr(i-1);
        end
        xCenterUm = (X_bl + W_arr / 2) / 1000;
        Coord_um_full = repmat(xCenterUm.', 1, Ny);     % Nx x Ny, varies along rows (X)
    else
        Y_bl = zeros(1, Ny);
        Y_bl(1) = -H_arr(1);
        for j = 2:Ny
            Y_bl(j) = Y_bl(j-1) + H_arr(j-1);
        end
        yCenterUm = (Y_bl + H_arr / 2) / 1000;
        Coord_um_full = repmat(yCenterUm, Nx, 1);       % Nx x Ny, varies along cols (Y)
    end

    roiCoordUm = Coord_um_full(monitor_indices);
    ref_um = mean(roiCoordUm);

    Target_Temp_full = baseTempK + gradientKPerUm * (Coord_um_full - ref_um);
    Target_Temp_ROI = Target_Temp_full(monitor_indices);
end

function Step1b_Extract_G_Matrices_6G_Decoupled()
%STEP1B_EXTRACT_G_MATRICES_6G_DECOUPLED Build G matrices from the 4-group full sweep.
%
% Step1a now solves all 4 symmetric control groups directly (edge_interior,
% edge_side, middle_interior, middle_side). This step interpolates the
% probe coordinates for each sweep case and stores both center- and base-
% offset G matrices.

import com.comsol.model.*
import com.comsol.model.util.*

powerWidthSettings = get_power_width_settings();
T_amb = powerWidthSettings.Ambient_Temp_K;
Nx_core = 42;    % Number of core columns
Ny_core = 42;    % Number of core rows
W_core = 1000;   % Core-cell width (nm)
H_core = 1000;   % Core-cell height (nm)
W_edge = 4000;
H_edge = 4000;

Nx = Nx_core + 2;
Ny = Ny_core + 2;

W_arr = ones(1, Nx) * W_core;
W_arr(1) = W_edge;
W_arr(end) = W_edge;

H_arr = ones(1, Ny) * H_core;
H_arr(1) = H_edge;
H_arr(end) = H_edge;

Base_P = 0.001; % W/um^2 probe density
Num_Pixels = Nx * Ny;
[masks, group_info, monitor_mask] = build_yrow_group_masks(Nx, Ny);
Num_Groups = numel(masks);
monitor_indices = find(monitor_mask);

targetSettings = get_thermal_target_settings();
Target_Times = targetSettings.Time_Step_s:targetSettings.Time_Step_s:targetSettings.Target_Time_s;
Num_Times = numel(Target_Times);

ModelFile = sprintf('ICM_Hetero_YRow_%dx%d.mph', Nx, Ny);
OutputFile = sprintf('G_Matrix_Hetero_YRow_%dx%d.mat', Nx, Ny);

fprintf('Loading base model: %s\n', ModelFile);
if ~exist(ModelFile, 'file')
    error('Cannot find %s. Run Step1a first.', ModelFile);
end
model = mphload(ModelFile);
fprintf('Base model loaded. Sweep covers %d control groups.\n', Num_Groups);

[X_centers, Y_centers] = build_pixel_centers(W_arr, H_arr);
Coords = build_probe_coords(X_centers, Y_centers, W_arr);

fprintf('Interpolating actual array coordinates: %d points, %d time samples.\n', ...
    size(Coords, 2), Num_Times);
tic;
Raw_Data = mphinterp(model, 'T', 'coord', Coords, ...
    'dataset', 'dset2', 't', Target_Times, 'outersolnum', 'all');
fprintf('Interpolation elapsed time: %.2f s.\n', toc);

[G_All, Num_Source_Groups] = reshape_probe_data(Raw_Data, T_amb, Base_P, ...
    Num_Times, 3, Num_Pixels);

if Num_Source_Groups ~= Num_Groups
    error('Expected %d sweep cases, but mphinterp returned %d solved group(s).', ...
        Num_Groups, Num_Source_Groups);
end

G_Data_Series = struct();
for tIdx = 1:Num_Times
    G_Data_Series(tIdx).time = Target_Times(tIdx);

    G_c_full = reshape(G_All(1, :, tIdx, :), Num_Pixels, Num_Groups);
    G_b_full = reshape(mean(G_All(2:3, :, tIdx, :), 1), Num_Pixels, Num_Groups);

    G_Data_Series(tIdx).G_center = G_c_full(monitor_indices, :);
    G_Data_Series(tIdx).G_base = G_b_full(monitor_indices, :);
    G_Data_Series(tIdx).G_center_full = G_c_full;
    G_Data_Series(tIdx).G_base_full = G_b_full;
end

save(OutputFile, 'G_Data_Series', 'masks', 'group_info', ...
    'monitor_indices', 'T_amb', 'W_arr', 'H_arr');
fprintf('Saved G matrix data to %s\n', OutputFile);
end

function [X_centers, Y_centers] = build_pixel_centers(W_arr, H_arr)
    Nx = numel(W_arr);
    Ny = numel(H_arr);

    X_bl = zeros(1, Nx);
    X_bl(1) = -W_arr(1);
    for i = 2:Nx
        X_bl(i) = X_bl(i-1) + W_arr(i-1);
    end
    X_centers = X_bl + W_arr / 2;

    Y_bl = zeros(1, Ny);
    Y_bl(1) = -H_arr(1);
    for j = 2:Ny
        Y_bl(j) = Y_bl(j-1) + H_arr(j-1);
    end
    Y_centers = Y_bl + H_arr / 2;
end

function Coords = build_probe_coords(X_centers, Y_values, W_arr)
    Nx = numel(X_centers);
    NyProbe = numel(Y_values);
    Coords = zeros(3, Nx * NyProbe * 3);

    idx = 0;
    for yIdx = 1:NyProbe
        for i = 1:Nx
            offsets_x = [0, -W_arr(i) / 2, W_arr(i) / 2];
            for k = 1:3
                idx = idx + 1;
                Coords(:, idx) = [X_centers(i) + offsets_x(k); ...
                    Y_values(yIdx); 110];
            end
        end
    end
end

function [G_Final, Num_Source_Groups] = reshape_probe_data( ...
    Raw_Data, T_amb, Base_P, Num_Times, Num_Probe_Positions, Num_Probe_Pixels)

    rawVector = Raw_Data(:);
    samplesPerGroup = Num_Times * Num_Probe_Positions * Num_Probe_Pixels;
    if mod(numel(rawVector), samplesPerGroup) ~= 0
        error(['Unexpected mphinterp output size. Got %d values, but each ' ...
            'solved group should provide %d values.'], ...
            numel(rawVector), samplesPerGroup);
    end

    Num_Source_Groups = numel(rawVector) / samplesPerGroup;
    G_4D = reshape((rawVector - T_amb) / Base_P, ...
        Num_Times, Num_Probe_Positions, Num_Probe_Pixels, Num_Source_Groups);
    G_Final = permute(G_4D, [2, 3, 1, 4]);
end

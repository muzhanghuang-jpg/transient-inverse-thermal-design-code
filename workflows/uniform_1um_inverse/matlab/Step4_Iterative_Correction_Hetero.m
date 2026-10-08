function [RMSE_ROI_meas, Max_Pos_Err_meas, Max_Neg_Err_meas] = ...
    Step4_Iterative_Correction_Hetero(iter, Gain)
%STEP4_ITERATIVE_CORRECTION_HETERO Closed-loop correction for Y-gradient target.

import com.comsol.model.*
import com.comsol.model.util.*

targetSettings = get_thermal_target_settings();
Target_Time = targetSettings.Target_Time_s;
Target_Temp = targetSettings.Base_Temp_K;
Target_Y_Gradient_K_per_um = targetSettings.Y_Gradient_K_per_um;
powerWidthSettings = get_power_width_settings();

Max_Contrast = powerWidthSettings.Max_Contrast;
rho_elec = powerWidthSettings.Rho_Elec_Ohm_m;
thickness = powerWidthSettings.Device_Thickness_m;
geometrySettings = get_geometry_constraint_settings();
Min_Wire_Width_nm = geometrySettings.Min_Wire_Width_nm;
Min_Side_Gap_nm = geometrySettings.Min_Side_Gap_nm;

fprintf('Loading base heterogeneous G matrix data...\n');
matFiles_G = dir('G_Matrix_Hetero_*.mat');
if isempty(matFiles_G)
    error('Cannot find G_Matrix_Hetero_*.mat in the current run directory.');
end
[~, latest_idx] = max([matFiles_G.datenum]);
G_Data_File = matFiles_G(latest_idx).name;
load(G_Data_File, 'G_Data_Series', 'masks', 'monitor_indices', ...
    'W_arr', 'H_arr');

Nx = length(W_arr);
Ny = length(H_arr);
Num_Pixels = Nx * Ny;

[Target_Temp_ROI, Target_Temp_full, Y_um_full, Target_Y_Reference_um] = ...
    build_gradient_target_profile(W_arr, H_arr, monitor_indices, ...
    Target_Temp, Target_Y_Gradient_K_per_um);

VerifyModelFile = sprintf('Verify_Model_Hetero_%dx%d_diode_%d.mph', Nx, Ny, iter);
OutFileName = sprintf('Correction_Step_Latest_Hetero_%dx%d_diode_%d.mat', Nx, Ny, iter + 1);
graphfilename = sprintf('Correction_Step_Latest_Hetero_%dx%d_diode_%d.pdf', Nx, Ny, iter + 1);

time_vec = [G_Data_Series.time];
[~, time_idx] = min(abs(time_vec - Target_Time));
G_c = G_Data_Series(time_idx).G_center;
G_b = G_Data_Series(time_idx).G_base;
Num_Groups = size(G_c, 2);

if iter == 0
    Latest_Correction_File = sprintf('Initial_Optimization_Solution_Hetero_%dx%d.mat', Nx, Ny);
    load(Latest_Correction_File, 'P_optimal_mW', 'Voltage_Use');
    P_Prev_W = P_optimal_mW / 1000;
else
    Latest_Correction_File = sprintf('Correction_Step_Latest_Hetero_%dx%d_diode_%d.mat', Nx, Ny, iter);
    load(Latest_Correction_File, 'P_New_mW', 'Voltage_Use');
    P_Prev_W = P_New_mW / 1000;
end

fprintf('Loading verification model: %s\n', VerifyModelFile);
if ~exist(VerifyModelFile, 'file')
    error('Cannot find verification model file: %s', VerifyModelFile);
end
try
    ModelUtil.remove('Model_Verify');
catch
end
% Always load the requested iter file. Rollback may evaluate iter N-1 while
% the in-memory Model_Verify still points to the rejected iter N model.
model = mphload(VerifyModelFile, 'Model_Verify');

X_bl = zeros(1, Nx);
X_bl(1) = -W_arr(1);
for i = 2:Nx
    X_bl(i) = X_bl(i-1) + W_arr(i-1);
end

Y_bl = zeros(1, Ny);
Y_bl(1) = -H_arr(1);
for j = 2:Ny
    Y_bl(j) = Y_bl(j-1) + H_arr(j-1);
end

Coords = zeros(3, Num_Pixels * 3);
idx = 0;
for j = 1:Ny
    for i = 1:Nx
        xc = X_bl(i) + W_arr(i) / 2;
        yc = Y_bl(j) + H_arr(j) / 2;
        offsets_x = [0, -W_arr(i) / 2, W_arr(i) / 2];
        for k = 1:3
            idx = idx + 1;
            Coords(:, idx) = [xc + offsets_x(k); yc; 110];
        end
    end
end

fprintf('Extracting measured temperature field at t = %.2f us...\n', Target_Time * 1e6);
tic;
Raw_Data = mphinterp(model, 'T', 'coord', Coords, 't', Target_Time);
fprintf('Extraction elapsed time: %.2f s.\n', toc);

T_measured_Raw = reshape(Raw_Data, 3, Num_Pixels);
T_b_measured_full = mean(T_measured_Raw(2:3, :), 1)';

T_b_measured = T_b_measured_full(monitor_indices);
T_map_measured = reshape(T_b_measured_full, Nx, Ny);
Error_map_measured = T_map_measured - Target_Temp_full;

Error_ROI_meas = T_b_measured - Target_Temp_ROI;
RMSE_ROI_meas = sqrt(mean(Error_ROI_meas .^ 2));
Max_Pos_Err_meas = max(Error_ROI_meas);
Max_Neg_Err_meas = min(Error_ROI_meas);

fprintf('\n--- Measured convergence against Y-gradient target ---\n');
fprintf('ROI measured RMSE: %.4f K\n', RMSE_ROI_meas);
fprintf('ROI max positive error: %+.4f K\n', Max_Pos_Err_meas);
fprintf('ROI max negative error: %+.4f K\n', Max_Neg_Err_meas);

fprintf('Computing power correction with Gain = %.3f...\n', Gain);
Target_Error = Error_ROI_meas;

C1 = G_b;
d1 = -Gain * Target_Error;
C2 = 0 * (G_c - G_b);
d2 = zeros(size(C2, 1), 1);
C_final = [C1; C2];
d_final = [d1; d2];

Constraints_List = build_yrow_group_adjacency_constraints();
num_constr = size(Constraints_List, 1);
A = zeros(num_constr, Num_Groups);
b = zeros(num_constr, 1);
for i = 1:num_constr
    g1 = Constraints_List(i, 1);
    g2 = Constraints_List(i, 2);
    A(i, g1) = 1;
    A(i, g2) = -Max_Contrast;
    b(i) = Max_Contrast * P_Prev_W(g2) - P_Prev_W(g1);
end

Aeq = [];
beq = [];
lb = -P_Prev_W;
ub = [];

options = optimoptions('lsqlin', 'Display', 'final');
dP_W = lsqlin(C_final, d_final, A, b, Aeq, beq, lb, ub, [], options);

P_New_W = P_Prev_W + dP_W;
P_New_mW = P_New_W * 1000;

power_density_Watts_New = zeros(Nx, Ny);
power_density_Watts_Prev = zeros(Nx, Ny);
for k = 1:Num_Groups
    power_density_Watts_New(masks{k}) = P_New_W(k);
    power_density_Watts_Prev(masks{k}) = P_Prev_W(k);
end

width_map_New_raw = zeros(Nx, Ny);
P_true_Watts_New = zeros(Nx, Ny);
P_true_Watts_Prev = zeros(Nx, Ny);

for i = 1:Nx
    for j = 1:Ny
        Area_um2 = (W_arr(i) / 1000) * (H_arr(j) / 1000);
        P_true_Watts_New(i, j) = power_density_Watts_New(i, j) * Area_um2;
        P_true_Watts_Prev(i, j) = power_density_Watts_Prev(i, j) * Area_um2;
    end
end

for i = 1:Nx
    row_power_total = sum(P_true_Watts_New(i, :));
    if row_power_total > 1e-12
        I_row = row_power_total / Voltage_Use;
        for j = 1:Ny
            P_pixel = P_true_Watts_New(i, j);
            if P_pixel > 1e-12
                R_pixel = P_pixel / (I_row ^ 2);
                L_meter = H_arr(j) * 1e-9;
                width_map_New_raw(i, j) = (rho_elec * L_meter) / (R_pixel * thickness) * 1e9;
            else
                width_map_New_raw(i, j) = 0;
            end
        end
    end
end

[Voltage_Use, width_map_New] = adapt_voltage_for_side_gap(width_map_New_raw, ...
    W_arr, Voltage_Use, Min_Wire_Width_nm, Min_Side_Gap_nm, ...
    sprintf('Correction width generation for iter %d', iter + 1));

save(OutFileName, 'P_New_mW', 'width_map_New', 'Voltage_Use', ...
    'power_density_Watts_New', 'P_true_Watts_New', 'Target_Time', ...
    'dP_W', 'Target_Temp', 'Target_Y_Gradient_K_per_um', ...
    'Target_Y_Reference_um', 'Target_Temp_ROI', 'Target_Temp_full', ...
    'Y_um_full');
fprintf('Saved correction data to: %s\n', OutFileName);

figure('Position', [100, 100, 1500, 950], 'Color', 'w');

subplot(2, 3, 1);
plot_smart_grid(Target_Temp_full, 'Target Temp (K)', '%.0f', 'jet');

subplot(2, 3, 2);
plot_smart_grid(T_map_measured, 'Measured Temp (K)', '%.0f', 'jet');
if ~all(isnan(T_map_measured(:)))
    clim([min(T_map_measured(:)), max(T_map_measured(:))]);
end

subplot(2, 3, 3);
plot_smart_grid(Error_map_measured, 'Measured Error (K)', '%+.1f', 'parula');

subplot(2, 3, 4);
dP_true_map_mW = (P_true_Watts_New - P_true_Watts_Prev) * 1000;
plot_smart_grid(dP_true_map_mW, 'True Power Correction (mW)', '%+.2f', 'parula');

subplot(2, 3, 5);
plot_smart_grid(P_true_Watts_New * 1000, 'New True Pixel Power (mW)', '%.2f', 'jet');
clim([0, max(P_true_Watts_New(:) * 1000)]);

subplot(2, 3, 6);
plot_smart_grid(width_map_New, sprintf('New Width @ %.0fV (nm)', Voltage_Use), '%.0f', 'parula');
clim([50, max(W_arr)]);

try
    exportgraphics(gcf, graphfilename, 'ContentType', 'vector');
    fprintf('Saved correction figure to: %s\n', graphfilename);
catch ME
    fprintf('\nNon-fatal warning: could not save PDF file: %s\n', graphfilename);
    fprintf('Reason: %s\n\n', ME.message);
end
end

function plot_smart_grid(data, title_str, num_fmt, cmap_name)
    im = imagesc(data');
    set(im, 'AlphaData', ~isnan(data'));
    axis equal tight;
    set(gca, 'YDir', 'normal');
    title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('X index (i)');
    ylabel('Y index (j)');

    try
        colormap(gca, cmap_name);
    catch
        colormap(gca, 'jet');
    end
    colorbar;
    hold on;

    [rows, cols] = size(data);
    for x = 0.5:1:rows+0.5
        line([x x], [0.5 cols+0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.1);
    end
    for y = 0.5:1:cols+0.5
        line([0.5 rows+0.5], [y y], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.1);
    end

    for r = 1:rows
        for c = 1:cols
            val = data(r, c);
            if isnan(val)
                continue;
            end
            if contains(title_str, 'Correction') && abs(val) < 1e-3
                continue;
            end
            if contains(title_str, 'Power') && abs(val) < 1e-3
                continue;
            end
            text(r, c, sprintf(num_fmt, val), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 2, 'Color', 'k');
        end
    end
end

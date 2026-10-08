function Step2b_optimization_Physical_WidthGen_Hetero()
%STEP2B_OPTIMIZATION_PHYSICAL_WIDTHGEN_HETERO Initial Y-gradient optimization.

targetSettings = get_thermal_target_settings();
Target_Time = targetSettings.Target_Time_s;
Target_Temp = targetSettings.Base_Temp_K;
Target_Y_Gradient_K_per_um = targetSettings.Y_Gradient_K_per_um;
GradientDirection = targetSettings.GradientDirection;
powerWidthSettings = get_power_width_settings();
T_ambient = powerWidthSettings.Ambient_Temp_K;
Voltage_Use = powerWidthSettings.Initial_Voltage_V;

Max_Contrast = powerWidthSettings.Max_Contrast;
rho_elec = powerWidthSettings.Rho_Elec_Ohm_m;
thickness = powerWidthSettings.Device_Thickness_m;
geometrySettings = get_geometry_constraint_settings();
Min_Wire_Width_nm = geometrySettings.Min_Wire_Width_nm;
Min_Side_Gap_nm = geometrySettings.Min_Side_Gap_nm;

matFiles = dir('G_Matrix_Hetero_*.mat');
if isempty(matFiles)
    error('Cannot find G_Matrix_Hetero_*.mat in the current run directory.');
end
[~, latest_idx] = max([matFiles.datenum]);
G_Data_File = matFiles(latest_idx).name;
fprintf('Loading G matrix data: %s\n', G_Data_File);
load(G_Data_File, 'G_Data_Series', 'masks', 'monitor_indices', 'W_arr', 'H_arr');

Nx = length(W_arr);
Ny = length(H_arr);
fprintf('Detected array size: Nx = %d, Ny = %d\n', Nx, Ny);

[Target_Temp_ROI, Target_Temp_full, Y_um_full, Target_Y_Reference_um] = ...
    build_gradient_target_profile(W_arr, H_arr, monitor_indices, ...
    Target_Temp, Target_Y_Gradient_K_per_um, GradientDirection);
fprintf('Target profile: %.2f K at ROI center, %+g K/um along +%s.\n', ...
    Target_Temp, Target_Y_Gradient_K_per_um, GradientDirection);

time_vec = [G_Data_Series.time];
[~, time_idx] = min(abs(time_vec - Target_Time));
fprintf('Optimizing target time: t = %.2f us\n', time_vec(time_idx) * 1e6);

G_c_calc = G_Data_Series(time_idx).G_center;
G_b_calc = G_Data_Series(time_idx).G_base;
Num_Groups = size(G_c_calc, 2);
fprintf('Detected control group count: %d\n', Num_Groups);

% Y-row block coalescing reduces wire-width steps in Y mode. In X mode the
% per-X-column grouping already enforces one width per column, so block
% coalescing is unnecessary and the expansion matrix becomes identity.
if strcmpi(GradientDirection, 'Y')
    blockSettings = get_yrow_block_settings();
    Y_Row_Block_Size = blockSettings.Y_Row_Block_Size;
    [group_expansion_matrix, yrow_block_info, row_to_yrow_block] = ...
        build_yrow_block_expansion(Ny, Y_Row_Block_Size);
    Num_Block_Groups = size(group_expansion_matrix, 2);
    if size(group_expansion_matrix, 1) ~= Num_Groups
        error('Y-row block expansion has %d rows, but G matrix has %d groups.', ...
            size(group_expansion_matrix, 1), Num_Groups);
    end
    fprintf('Y_Row_Block_Size = %d; optimizing %d block control groups.\n', ...
        Y_Row_Block_Size, Num_Block_Groups);
else
    Y_Row_Block_Size = 0;
    yrow_block_info = struct();
    row_to_yrow_block = [];
    group_expansion_matrix = eye(Num_Groups);
    Num_Block_Groups = Num_Groups;
    fprintf('X-gradient mode: per-X-column grouping; block coalescing skipped.\n');
end

C1 = G_b_calc;
d1 = Target_Temp_ROI - T_ambient;
C2 = 0 * (G_c_calc - G_b_calc);
d2 = zeros(size(C2, 1), 1);
C_final = [C1; C2];
d_final = [d1; d2];
C_solve = C_final * group_expansion_matrix;

Nprimary = Num_Groups / 2;   % primary-axis block count (Ny in Y mode, Nx in X mode)
Constraints_List = build_yrow_group_adjacency_constraints(Nprimary);
num_constr = size(Constraints_List, 1);
A_full = zeros(num_constr, Num_Groups);
b = zeros(num_constr, 1);
for i = 1:num_constr
    A_full(i, Constraints_List(i, 1)) = 1;
    A_full(i, Constraints_List(i, 2)) = -Max_Contrast;
end
A = A_full * group_expansion_matrix;

Aeq = [];
beq = [];
lb = zeros(Num_Block_Groups, 1);

if strcmpi(GradientDirection, 'X') && powerWidthSettings.X_Smoothing_Lambda > 0
    lambda_smooth = powerWidthSettings.X_Smoothing_Lambda;
    D_smooth = build_column_smoothness_matrix(Nprimary);
    C_solve = [C_solve; lambda_smooth * (D_smooth * group_expansion_matrix)];
    d_final = [d_final; zeros(size(D_smooth, 1), 1)];
    fprintf('Applied X-column smoothness regularisation: lambda = %.4f, %d penalty rows.\n', ...
        lambda_smooth, size(D_smooth, 1));
end

fprintf('Solving initial lsqlin problem for %s-gradient control groups...\n', GradientDirection);
options = optimoptions('lsqlin', 'Display', 'final-detailed');
P_block_scale = lsqlin(C_solve, d_final, A, b, Aeq, beq, lb, [], [], options);
P_scale = group_expansion_matrix * P_block_scale;

P_absolute_Watts = P_scale;
P_optimal_mW = P_absolute_Watts * 1000;
P_block_mW = P_block_scale * 1000;

T_ROI_pred = T_ambient + G_b_calc * P_scale;
Error_ROI_pred = T_ROI_pred - Target_Temp_ROI;
RMSE_ROI_pred = sqrt(mean(Error_ROI_pred .^ 2));
Max_Pos_Err_pred = max(Error_ROI_pred);
Max_Neg_Err_pred = min(Error_ROI_pred);

fprintf('\n--- Predicted convergence against %s-gradient target ---\n', GradientDirection);
fprintf('ROI predicted RMSE: %.2f K\n', RMSE_ROI_pred);
fprintf('ROI max positive error: %+.2f K\n', Max_Pos_Err_pred);
fprintf('ROI max negative error: %+.2f K\n\n', Max_Neg_Err_pred);

power_density_Watts = zeros(Nx, Ny);
for k = 1:Num_Groups
    power_density_Watts(masks{k}) = P_absolute_Watts(k);
end

width_map_raw = zeros(Nx, Ny);
P_true_Watts = zeros(Nx, Ny);

for i = 1:Nx
    for j = 1:Ny
        W_um = W_arr(i) / 1000;
        H_um = H_arr(j) / 1000;
        P_true_Watts(i, j) = power_density_Watts(i, j) * W_um * H_um;
    end
end

for i = 1:Nx
    row_power_total = sum(P_true_Watts(i, :));
    if row_power_total > 1e-12
        I_row = row_power_total / Voltage_Use;
        for j = 1:Ny
            P_pixel = P_true_Watts(i, j);
            if P_pixel > 1e-12
                R_pixel = P_pixel / (I_row ^ 2);
                L_meter = H_arr(j) * 1e-9;
                width_map_raw(i, j) = (rho_elec * L_meter) / (R_pixel * thickness) * 1e9;
            else
                width_map_raw(i, j) = 0;
            end
        end
    end
end

[Voltage_Use, width_map_final] = adapt_voltage_for_side_gap(width_map_raw, ...
    W_arr, Voltage_Use, Min_Wire_Width_nm, Min_Side_Gap_nm, ...
    'Initial optimization width generation');

if isfield(G_Data_Series, 'G_base_full')
    G_b_full = G_Data_Series(time_idx).G_base_full;
    T_map_full = reshape(T_ambient + G_b_full * P_scale, Nx, Ny);
else
    T_map_full = NaN(Nx, Ny);
end
Error_map_full = T_map_full - Target_Temp_full;

OutFileName = sprintf('Initial_Optimization_Solution_Hetero_%dx%d.mat', Nx, Ny);
save(OutFileName, 'P_optimal_mW', 'width_map_final', 'Voltage_Use', ...
    'power_density_Watts', 'P_true_Watts', 'Target_Time', 'A', 'b', ...
    'P_scale', 'P_block_scale', 'P_block_mW', 'Y_Row_Block_Size', ...
    'yrow_block_info', 'row_to_yrow_block', 'group_expansion_matrix', ...
    'Target_Temp', 'Target_Y_Gradient_K_per_um', 'GradientDirection', ...
    'Target_Y_Reference_um', 'Target_Temp_ROI', 'Target_Temp_full', ...
    'Y_um_full');
fprintf('Saved initial optimization data to: %s\n', OutFileName);

figure('Position', [50, 50, 1600, 900], 'Color', 'w');
subplot(2, 3, 1);
plot_smart_grid(P_true_Watts * 1000, 'True Pixel Power (mW)', '%.2f', 'jet');
subplot(2, 3, 2);
plot_smart_grid(width_map_final, sprintf('Width @ %.0fV (nm)', Voltage_Use), '%.0f', 'parula');
subplot(2, 3, 3);
plot_smart_grid(Target_Temp_full, 'Target Temp (K)', '%.0f', 'jet');
subplot(2, 3, 4);
plot_smart_grid(T_map_full, sprintf('Predicted Temp (K) @ %.0fus', Target_Time * 1e6), '%.0f', 'jet');
subplot(2, 3, 5);
plot_smart_grid(Error_map_full, 'Predicted Error (K)', '%+.1f', 'parula');
subplot(2, 3, 6);
plot_smart_grid(reshape_vector_to_map(Error_ROI_pred, monitor_indices, Nx, Ny), ...
    'ROI Error (K)', '%+.1f', 'parula');

filename = sprintf('Optimization_Result_Hetero_%dx%d.pdf', Nx, Ny);
exportgraphics(gcf, filename, 'ContentType', 'vector');
fprintf('Saved initial optimization figure to: %s\n', filename);
end

function out = reshape_vector_to_map(values, indices, Nx, Ny)
    out = NaN(Nx, Ny);
    out(indices) = values;
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
            if contains(title_str, 'Power') && abs(val) < 1e-4
                continue;
            end
            text(r, c, sprintf(num_fmt, val), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 2, 'Color', 'k');
        end
    end
end

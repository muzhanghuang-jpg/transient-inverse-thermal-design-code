function Step2b_optimization_Physical_WidthGen_Hetero
% =========================================================================
% Step 2b: optimize group power densities and convert them to filament widths.
% 
% Six-group rectangular design with separate top and bottom groups.
% Read group masks and cell dimensions from the response-matrix file.
% Top and bottom powers are not constrained to be equal.
% Display the y-axis in increasing physical-coordinate order.
% =========================================================================

% --- 1. Parameters ---
Target_Time = 13E-6;  % Target time (s)
Target_Temp = 900;    % Target temperature (K)
T_ambient   = 293.15; % Ambient temperature (K)
Voltage_Use = 30;     % Applied voltage (V)
 
Max_Contrast = 4;     % Upper bound for the constrained group power-density ratios
rho_elec = 1 / (75272/6); % Resistivity used for the initial width conversion (ohm m)
thickness = 220e-9;       % Silicon thickness (m)
 
% --- 2. Load the response data ---
matFiles = dir('G_Matrix_Hetero_*.mat');
if isempty(matFiles)
    error('Response data are missing; expected G_Matrix_Hetero_*.mat in the current directory.');
else
    [~, latest_idx] = max([matFiles.datenum]);
    G_Data_File = matFiles(latest_idx).name; 
    fprintf('Load the latest response data: %s ...\n', G_Data_File);
end
load(G_Data_File); % Includes G_Data_Series, masks, W_arr, and H_arr.

Nx = length(W_arr);
Ny = length(H_arr);
fprintf('Grid dimensions: Nx = %d, Ny = %d\n', Nx, Ny);

% Select the response time closest to the target time.
time_vec = [G_Data_Series.time];
[~, time_idx] = min(abs(time_vec - Target_Time));
fprintf('Optimize the sampled response at t = %.2f us.\n', time_vec(time_idx)*1e6);

% Target-region response matrices; each column corresponds to a virtual group.
G_c_calc = G_Data_Series(time_idx).G_center; 
G_b_calc = G_Data_Series(time_idx).G_base;   
Num_Groups = size(G_c_calc, 2); 
fprintf('Number of virtual groups: %d\n', Num_Groups);

% --- 3. Least-squares objective ---
Target_Rise = Target_Temp - T_ambient;
C1 = G_b_calc;
d1 = ones(size(G_b_calc, 1), 1) * Target_Rise;
C2 = 0 * (G_c_calc - G_b_calc); % The center-to-boundary variation penalty is disabled.
d2 = zeros(size(C2, 1), 1);
C_final = [C1; C2];
d_final = [d1; d2];

% --- 4. Constraints ---
% 4.1 Group power-density ratio constraints (A*x <= b)
Constraints_List = build_constraints_list(masks, Nx, Ny, Num_Groups);
num_constr = size(Constraints_List, 1);
A = zeros(num_constr, Num_Groups); b = zeros(num_constr, 1);
for i = 1:num_constr
    A(i, Constraints_List(i, 1)) = 1;
    A(i, Constraints_List(i, 2)) = -Max_Contrast;
end

% 4.2 Equality constraints (Aeq*x = beq)
% Leave these empty so the top and bottom groups can vary independently.
Aeq = []; beq = [];

% --- 5. Constrained solve ---
fprintf('Solve the six-group constrained least-squares problem.\n');
options = optimoptions('lsqlin', 'Display', 'final-detailed');
lb = zeros(Num_Groups, 1);
P_scale = lsqlin(C_final, d_final, A, b, Aeq, beq, lb, [], [], options);

P_absolute_Watts = P_scale; 
P_optimal_mW = P_absolute_Watts * 1000;

% --- Evaluate the linear prediction at the optimization points ---
T_ROI_pred = T_ambient + G_b_calc * P_scale; 
Error_ROI_pred = T_ROI_pred - Target_Temp;
RMSE_ROI_pred = sqrt(mean(Error_ROI_pred.^2));
Max_Pos_Err_pred = max(Error_ROI_pred); 
Max_Neg_Err_pred = min(Error_ROI_pred); 

fprintf('\n--- Sampled errors predicted by the response matrix ---\n');
fprintf('Predicted target-region RMSE: %.2f K\n', RMSE_ROI_pred);
fprintf('Maximum positive sampled error: %+.2f K\n', Max_Pos_Err_pred);
fprintf('Maximum negative sampled error: %+.2f K\n\n', Max_Neg_Err_pred);

% --- 6. Cell powers and filament widths ---
power_density_Watts = zeros(Nx, Ny);
for k = 1:Num_Groups
    power_density_Watts(masks{k}) = P_absolute_Watts(k);
end

width_map_final = zeros(Nx, Ny); 
P_true_Watts = zeros(Nx, Ny); 

for i = 1:Nx
    for j = 1:Ny
        W_um = W_arr(i) / 1000; 
        H_um = H_arr(j) / 1000;
        Area_um2 = W_um * H_um;
        % Cell power equals group power density multiplied by physical cell area.
        P_true_Watts(i, j) = power_density_Watts(i, j) * Area_um2;
    end
end

for i = 1:Nx
    row_power_total = sum(P_true_Watts(i, :)); 
    if row_power_total > 1e-12
        I_row = row_power_total / Voltage_Use; 
        for j = 1:Ny
            P_pixel = P_true_Watts(i, j);
            if P_pixel > 1e-12
                R_pixel = P_pixel / (I_row^2);
                L_meter = H_arr(j) * 1e-9; 
                W_calc = (rho_elec * L_meter) / (R_pixel * thickness);
                
                W_max_meter = W_arr(i) * 1e-9;
                W_clamped = min(max(W_calc, 50e-9), W_max_meter); 
                width_map_final(i, j) = W_clamped * 1e9; 
            else
                 width_map_final(i, j) = 0; 
            end
        end
    end
end

% 6.3 Reconstruct sampled temperatures across the full grid.
if isfield(G_Data_Series, 'G_base_full')
    G_b_full = G_Data_Series(time_idx).G_base_full; 
    T_map_full = reshape(T_ambient + G_b_full * P_scale, Nx, Ny);
else
    T_map_full = NaN(Nx, Ny);
end

% --- 7. Save ---
OutFileName = sprintf('Initial_Optimization_Solution_Hetero_%dx%d.mat', Nx, Ny);
save(OutFileName, 'P_optimal_mW', 'width_map_final', 'Voltage_Use', 'power_density_Watts', 'P_true_Watts', 'Target_Time', 'A', 'b', 'P_scale', 'Target_Temp');
fprintf('Initial optimization data saved: %s\n', OutFileName);

% --- 8. Plots ---
figure('Position', [50, 50, 1600, 500], 'Color', 'w');
subplot(1, 3, 1);
plot_smart_grid(P_true_Watts * 1000, 'True Pixel Power (mW)', '%.2f', 'jet');
subplot(1, 3, 2);
plot_smart_grid(width_map_final, sprintf('Width @ %.0fV (nm)', Voltage_Use), '%.0f', 'parula');
subplot(1, 3, 3);
plot_smart_grid(T_map_full, sprintf('Predicted Temp (K) @ %.0fus', Target_Time*1e6), '%.0f', 'jet');

filename = sprintf('Optimization_Result_Hetero_%dx%d.pdf', Nx, Ny);
exportgraphics(gcf, filename, 'ContentType', 'vector');
end

% --- Plotting helper ---
function plot_smart_grid(data, title_str, num_fmt, cmap_name)
    im = imagesc(data'); 
    set(im, 'AlphaData', ~isnan(data')); 
    axis equal tight;
    
    % Physical y-coordinates increase upward.
    set(gca, 'YDir', 'normal'); 
    
    title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Row Index (i)'); ylabel('Col Index (j)');
    try colormap(gca, cmap_name); catch; colormap(gca, 'jet'); end
    colorbar; hold on;
    
    [rows, cols] = size(data);
    for x = 0.5:1:rows+0.5, line([x x], [0.5 cols+0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.1); end
    for y = 0.5:1:cols+0.5, line([0.5 rows+0.5], [y y], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.1); end
    
    c_limits = clim; c_min = c_limits(1); c_range = c_limits(2) - c_min;
    if c_range == 0, c_range = 1; end
    n_colors = 256;
    
    for r = 1:rows
        for c = 1:cols
            val = data(r, c);
            if isnan(val) || (contains(title_str, 'Power') && val < 1e-4), continue; end
            text(r, c, sprintf(num_fmt, val), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 2, 'Color', 'k');
        end
    end
end

% --- Constraint-list helper ---
function Constraints_List = build_constraints_list(masks, Nx, Ny, Num_Groups)
    Constraints_List = [];
    for r = 1:Nx
        g_row = []; for k=1:Num_Groups, if any(masks{k}(r,:)), g_row=[g_row, k]; end; end
        g_row = unique(g_row);
        for i=1:length(g_row), for j=1:length(g_row), if i~=j, Constraints_List=[Constraints_List; g_row(i), g_row(j)]; end; end; end
    end
    for c = 1:Ny
        g_col = []; for k=1:Num_Groups, if any(masks{k}(:,c)), g_col=[g_col, k]; end; end
        g_col = unique(g_col);
        for i=1:length(g_col), for j=1:length(g_col), if i~=j, Constraints_List=[Constraints_List; g_col(i), g_col(j)]; end; end; end
    end
    Constraints_List = unique(Constraints_List, 'rows');
end

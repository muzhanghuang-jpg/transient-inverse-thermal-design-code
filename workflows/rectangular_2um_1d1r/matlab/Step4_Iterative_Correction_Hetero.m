function [RMSE_ROI_meas, Max_Pos_Err_meas, Max_Neg_Err_meas] = Step4_Iterative_Correction_Hetero(iter, Gain)
% =========================================================================
% Step 4: update group power densities using nonlinear FEM temperature errors.
% 
% Grid dimensions W_arr and H_arr are loaded from the response data.
% Extract samples in MATLAB column-major order: j outside, i inside.
% Convert group densities to cell powers before calculating filament widths.
% Read the initial solution or the preceding accepted correction.
% The correction uses separate top and bottom virtual groups.
% =========================================================================

import com.comsol.model.*
import com.comsol.model.util.*

% --- 1. Parameters ---
Target_Time = 13E-6; % Must match the response and verification target time.
Target_Temp = 900;                        
                             
Max_Contrast = 4;      

rho_elec = 1 / (75272/6); % Resistivity used for width conversion (ohm m)
thickness = 220e-9; 

% --- 2. Load response and previous power data ---
fprintf('Load the nonuniform-grid response matrices.\n');
matFiles_G = dir('G_Matrix_Hetero_*.mat');
if isempty(matFiles_G), error('Response data are missing; expected G_Matrix_Hetero_*.mat.'); end
G_Data_File = matFiles_G(1).name;
load(G_Data_File, 'G_Data_Series', 'masks', 'monitor_indices', 'T_amb', 'W_arr', 'H_arr'); 

Nx = length(W_arr);
Ny = length(H_arr);
Num_Pixels = Nx * Ny;

% Construct model and output filenames.

VerifyModelFile = sprintf('Verify_Model_Hetero_%dx%d_diode_%d.mph', Nx, Ny, iter);
OutFileName = sprintf('Correction_Step_Latest_Hetero_%dx%d_diode_%d.mat', Nx, Ny, iter+1);
graphfilename = sprintf('Correction_Step_Latest_Hetero_%dx%d_diode_%d.pdf', Nx, Ny, iter+1);

% Number of virtual groups.
time_vec = [G_Data_Series.time];
[~, time_idx] = min(abs(time_vec - Target_Time));
G_c = G_Data_Series(time_idx).G_center;
G_b = G_Data_Series(time_idx).G_base;

Num_Groups = size(G_c, 2);


% Read the preceding group power-density values.
if iter == 0
    Latest_Correction_File = sprintf('Initial_Optimization_Solution_Hetero_%dx%d.mat', Nx, Ny);
    load(Latest_Correction_File, 'P_optimal_mW', 'Voltage_Use');
    P_Prev_W = P_optimal_mW / 1000; % Convert the initial stored values.
else
    Latest_Correction_File = sprintf('Correction_Step_Latest_Hetero_%dx%d_diode_%d.mat', Nx, Ny, iter);
    load(Latest_Correction_File, 'P_New_mW', 'Voltage_Use');
    P_Prev_W = P_New_mW / 1000; % Convert the preceding corrected values.
end


% --- 3. Extract nonlinear FEM temperatures ---
fprintf('Load the verification model: %s ...\n', VerifyModelFile);
if ~exist(VerifyModelFile, 'file'), error('Verification model file is missing.'); end
model = mphload(VerifyModelFile);

% Construct sampling coordinates on the nonuniform grid.
X_bl = zeros(1, Nx); X_bl(1) = -W_arr(1); 
for i = 2:Nx, X_bl(i) = X_bl(i-1) + W_arr(i-1); end
Y_bl = zeros(1, Ny); Y_bl(1) = -H_arr(1); 
for j = 2:Ny, Y_bl(j) = Y_bl(j-1) + H_arr(j-1); end

Coords = zeros(3, Num_Pixels * 3);
idx = 0;

% Keep sample ordering consistent with the response-matrix extraction.
for j = 1:Ny
    for i = 1:Nx
        xc = X_bl(i) + W_arr(i)/2;
        yc = Y_bl(j) + H_arr(j)/2;
        offsets_x = [0, -W_arr(i)/2, W_arr(i)/2]; 
        for k = 1:3
            idx = idx + 1;
            Coords(:, idx) = [xc + offsets_x(k); yc; 110];
        end
    end
end

fprintf('Extract nonlinear FEM temperatures at t = %.2f us.\n', Target_Time*1e6);
tic;
Raw_Data = mphinterp(model, 'T', 'coord', Coords, 't', Target_Time);
fprintf('Extraction time: %.2f s.\n', toc);

% Separate center samples and the mean transverse-boundary temperatures.
T_measured_Raw = reshape(Raw_Data, 3, Num_Pixels);
T_c_measured_full = T_measured_Raw(1, :)';
T_b_measured_full = mean(T_measured_Raw(2:3, :), 1)';

% Retain target-region temperatures for optimization.
T_c_measured = T_c_measured_full(monitor_indices);
T_b_measured = T_b_measured_full(monitor_indices);

% Reshape sampled temperatures across the full grid for plotting.
T_map_measured = reshape(T_b_measured_full, Nx, Ny);

% --- Evaluate nonlinear sampled errors ---
Error_ROI_meas = T_b_measured - Target_Temp;
RMSE_ROI_meas = sqrt(mean(Error_ROI_meas.^2));

% Maximum positive and negative errors at the optimization points.
Max_Pos_Err_meas = max(Error_ROI_meas); 
Max_Neg_Err_meas = min(Error_ROI_meas); 



% --- 4. Constrained power-density correction ---
fprintf('Calculate the group power-density correction.\n');

% Error relative to the prescribed target at the transverse-boundary samples.
Target_Error = T_b_measured - Target_Temp;

% Set the predicted temperature increment to -Gain * Error.
C1 = G_b;
d1 = -Gain * Target_Error;

% The center-to-boundary variation penalty remains disabled.
C2 = 0 * (G_c - G_b);
d2 = zeros(size(C2, 1), 1);
C_final = [C1; C2];
d_final = [d1; d2];

% 4.1 Group power-density ratio constraints
Constraints_List = build_constraints_list(masks, Nx, Ny, Num_Groups);
num_constr = size(Constraints_List, 1);
A = zeros(num_constr, Num_Groups); b = zeros(num_constr, 1);
for i = 1:num_constr
    g1 = Constraints_List(i, 1); g2 = Constraints_List(i, 2);
    A(i, g1) = 1; A(i, g2) = -Max_Contrast;
    % Shift the incremental constraint bounds using the previous absolute densities.
    b(i) = Max_Contrast * P_Prev_W(g2) - P_Prev_W(g1);
end

% 4.2 Equality constraints
% Leave these empty so the top and bottom groups can vary independently.
Aeq = []; beq = [];

% 4.3 Bounds prevent corrected group powers from becoming negative.
lb = -P_Prev_W; 
ub = [];

options = optimoptions('lsqlin', 'Display', 'final');
% Pass empty equality-constraint matrices to lsqlin.
dP_W = lsqlin(C_final, d_final, A, b, Aeq, beq, lb, ub, [], options);

% --- 5. Corrected cell powers and filament widths ---
P_New_W = P_Prev_W + dP_W;
P_New_mW = P_New_W * 1000;

power_density_Watts_New = zeros(Nx, Ny);
power_density_Watts_Prev = zeros(Nx, Ny);
for k = 1:Num_Groups
    power_density_Watts_New(masks{k}) = P_New_W(k);
    power_density_Watts_Prev(masks{k}) = P_Prev_W(k);
end

width_map_New = zeros(Nx, Ny); 
P_true_Watts_New = zeros(Nx, Ny); 
P_true_Watts_Prev = zeros(Nx, Ny);

for i = 1:Nx
    for j = 1:Ny
        width_map_New(i, j) = W_arr(i) * 1e9; % Initialize to the physical cell width.
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
                R_pixel = P_pixel / (I_row^2);
                L_meter = H_arr(j) * 1e-9;
                W_calc = (rho_elec * L_meter) / (R_pixel * thickness);
                
                W_max_meter = W_arr(i) * 1e-9;
                W_clamped = min(max(W_calc, 50e-9), W_max_meter); 
                width_map_New(i, j) = W_clamped * 1e9; 
            else
                width_map_New(i, j) = 0; 
            end
        end
    end
end

% --- 6. Save the corrected state ---

save(OutFileName, 'P_New_mW', 'width_map_New', 'Voltage_Use', 'power_density_Watts_New', 'P_true_Watts_New', 'Target_Time', 'dP_W');
fprintf('Corrected linewidth data saved: %s\n', OutFileName);

% --- 7. Verification plots ---
figure('Position', [100, 100, 1200, 1000], 'Color', 'w');

% Nonlinear FEM temperatures at cell centers.
subplot(2, 2, 1);
plot_smart_grid(T_map_measured, 'Measured Temp (K) before Corr', '%.0f', 'jet');
if ~all(isnan(T_map_measured(:))), clim([min(T_map_measured(:)), max(T_map_measured(:))]); end

% Change in cell power.
dP_true_map_mW = (P_true_Watts_New - P_true_Watts_Prev) * 1000;
subplot(2, 2, 2);
plot_smart_grid(dP_true_map_mW, 'True Power Correction (mW)', '%+.2f', 'parula');

% Corrected cell powers.
subplot(2, 2, 3);
plot_smart_grid(P_true_Watts_New * 1000, 'New True Pixel Power (mW)', '%.2f', 'jet');
clim([0, max(P_true_Watts_New(:)*1000)]);

% Corrected filament widths.
subplot(2, 2, 4);
plot_smart_grid(width_map_New, sprintf('New Width @ %.0fV (nm)', Voltage_Use), '%.0f', 'parula');
clim([50, max(W_arr)]);


try
    exportgraphics(gcf, graphfilename, 'ContentType', 'vector');
    fprintf('Figure saved: %s\n', graphfilename);
catch ME
    fprintf('\nWarning: could not save PDF: %s\n', graphfilename);
    fprintf('The output file may be open in another application.\n');
    fprintf('Figure output was skipped; saved MAT data and correction states are unaffected.\n\n');
end

end

% --- Plotting and constraint helpers ---
function plot_smart_grid(data, title_str, num_fmt, cmap_name)
    im = imagesc(data'); 
    set(im, 'AlphaData', ~isnan(data')); 
    axis equal tight;
    set(gca, 'YDir', 'normal')
    title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Matrix Row Index'); ylabel('Matrix Column Index');
    
    try cmap = colormap(gca, cmap_name); catch, cmap = colormap(gca, 'jet'); end
    colorbar; hold on;
    
    [rows, cols] = size(data);
    for x = 0.5 : 1 : rows+0.5, line([x x], [0.5 cols+0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.1); end
    for y = 0.5 : 1 : cols+0.5, line([0.5 rows+0.5], [y y], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.1); end
    
    c_limits = clim; c_min = c_limits(1); c_range = c_limits(2) - c_min;
    if c_range == 0, c_range = 1; end
    n_colors = size(cmap, 1);
    
    for r = 1:rows
        for c = 1:cols
            val = data(r, c);
            if isnan(val), continue; end
            
            if contains(title_str, 'Correction') && abs(val) < 1e-3, continue; end
            if contains(title_str, 'Power') && abs(val) < 1e-3, continue; end
            
            idx = floor( (val - c_min) / c_range * (n_colors-1) ) + 1;
            idx = max(1, min(idx, n_colors));
            
            text(r, c, sprintf(num_fmt, val), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 2, 'Color', 'k');
        end
    end
end

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

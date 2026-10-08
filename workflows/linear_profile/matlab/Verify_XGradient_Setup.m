function Verify_XGradient_Setup()
%VERIFY_XGRADIENT_SETUP Dry-run sanity check for X-gradient configuration.
%
% Validates that the helper layer (settings, masks, adjacency, target
% profile) returns shapes and indices consistent with GradientDirection = 'X'.
% No COMSOL connection required; runs in seconds.

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

fprintf('\n==================================================\n');
fprintf('X-gradient dry-run sanity check\n');
fprintf('==================================================\n\n');

% --- 1. Settings -----------------------------------------------------------
settings = get_thermal_target_settings();
fprintf('[1] get_thermal_target_settings\n');
fprintf('    GradientDirection      = %s\n', settings.GradientDirection);
fprintf('    Y_Gradient_K_per_um    = %g K/um\n', settings.Y_Gradient_K_per_um);
fprintf('    Base_Temp_K            = %g K\n', settings.Base_Temp_K);
assert(strcmpi(settings.GradientDirection, 'X'), ...
    'GradientDirection should be X for v2 default');

% --- 2. Array dimensions (matching Step1a/Step1b) --------------------------
Nx_core = 31; Ny_core = 21;
Nx = Nx_core + 2;     % 33
Ny = Ny_core + 2;     % 23
W_core = 2000; H_core = 2000; W_edge = 4000; H_edge = 4000;
W_arr = ones(1, Nx) * W_core; W_arr(1) = W_edge; W_arr(end) = W_edge;
H_arr = ones(1, Ny) * H_core; H_arr(1) = H_edge; H_arr(end) = H_edge;

% --- 3. Group masks --------------------------------------------------------
[masks, group_info, monitor_mask] = build_yrow_group_masks(Nx, Ny, settings.GradientDirection);
Num_Groups = numel(masks);
fprintf('\n[2] build_yrow_group_masks (X mode)\n');
fprintf('    Num_Groups             = %d (expected 2*Nx = %d)\n', Num_Groups, 2 * Nx);
assert(Num_Groups == 2 * Nx);

fprintf('    Group 1 (col 1 int.)   covers %d pixels (expected Ny-2 = %d)\n', ...
    sum(masks{1}(:)), Ny - 2);
assert(sum(masks{1}(:)) == Ny - 2);
assert(all(masks{1}(1, 2:Ny-1)), 'group 1 should mark col=1, j=2..Ny-1');

fprintf('    Group 2 (col 1 TB-edge) covers %d pixels (expected 2)\n', sum(masks{2}(:)));
assert(sum(masks{2}(:)) == 2);
assert(masks{2}(1, 1) && masks{2}(1, Ny));

fprintf('    Group %d (col Nx int.)  covers %d pixels\n', ...
    2*Nx-1, sum(masks{2*Nx-1}(:)));
assert(all(masks{2*Nx-1}(Nx, 2:Ny-1)));

assert(strcmp(group_info(1).type, 'interior'));
assert(group_info(1).col_i == 1);
assert(strcmp(group_info(2).type, 'side_tb'));
assert(strcmp(group_info(2*Nx).type, 'side_tb'));

fprintf('    monitor_mask covers     %d pixels (expected (Nx-2)*(Ny-2) = %d)\n', ...
    sum(monitor_mask(:)), (Nx-2)*(Ny-2));
assert(sum(monitor_mask(:)) == (Nx-2)*(Ny-2));

% --- 4. Fast-mode case selection -------------------------------------------
Center_Index = ceil(Nx / 2);
Center_Group_Indices = [2*Center_Index - 1, 2*Center_Index];
Edge_And_Center = [1, 2, Center_Group_Indices, 2*Nx - 1, 2*Nx];
fprintf('\n[3] fast-mode case ids (X mode, Nx = %d)\n', Nx);
fprintf('    Center_Index           = %d\n', Center_Index);
fprintf('    Edge_And_Center        = [%s]\n', num2str(Edge_And_Center));
expected_ids = [1, 2, 33, 34, 65, 66];
assert(isequal(Edge_And_Center, expected_ids), ...
    sprintf('case ids mismatch: %s vs %s', mat2str(Edge_And_Center), mat2str(expected_ids)));

% --- 5. Adjacency constraints ----------------------------------------------
Nprimary = Num_Groups / 2;
constraints = build_yrow_group_adjacency_constraints(Nprimary);
fprintf('\n[4] build_yrow_group_adjacency_constraints(Nprimary = %d)\n', Nprimary);
fprintf('    constraints rows       = %d\n', size(constraints, 1));
expected_constr = 6 * (Nprimary - 1) + 2;
fprintf('    expected               = %d (6*(N-1) + 2)\n', expected_constr);
assert(size(constraints, 1) == expected_constr);

% --- 6. Target profile -----------------------------------------------------
monitor_indices = find(monitor_mask);
[T_ROI, T_full, coord_um, ref_um] = build_gradient_target_profile( ...
    W_arr, H_arr, monitor_indices, ...
    settings.Base_Temp_K, settings.Y_Gradient_K_per_um, settings.GradientDirection);
fprintf('\n[5] build_gradient_target_profile (X mode)\n');
fprintf('    T_full size            = %s\n', mat2str(size(T_full)));
fprintf('    coord_um size          = %s\n', mat2str(size(coord_um)));
fprintf('    ROI ref coord          = %.3f um\n', ref_um);

% In X mode, coord_um and T_full vary along X (rows) and are constant along Y (cols).
coord_std_along_y = std(coord_um, 0, 2);
T_std_along_y = std(T_full, 0, 2);
fprintf('    max std(coord_um) along Y = %.6g (expected 0 in X mode)\n', max(coord_std_along_y));
fprintf('    max std(T_full)   along Y = %.6g K (expected 0 in X mode)\n', max(T_std_along_y));
assert(max(coord_std_along_y) < 1e-9);
assert(max(T_std_along_y) < 1e-9);

T_min = min(T_full(:)); T_max = max(T_full(:));
fprintf('    T_full range           = [%.2f, %.2f] K (span %.2f K)\n', ...
    T_min, T_max, T_max - T_min);

% --- 7. Summary ------------------------------------------------------------
fprintf('\n==================================================\n');
fprintf('All %d helper checks PASSED. X-gradient layer is consistent.\n', 5);
fprintf('==================================================\n');
end

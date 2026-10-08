function Step1b_Extract_G_Matrices_6G_Decoupled
% =========================================================================
% Step 1b: extract the six-group transient thermal response matrices.
% Grid dimensions and group order must match Step1a.
% =========================================================================
import com.comsol.model.*
import com.comsol.model.util.*

% --- 1. Parameters matching the response model ---
T_amb = 293.15; 
Nx_core = 31; Ny_core = 21;
W_core = 2000; H_core = 2000;
W_edge = 4000; H_edge = 4000;

Nx = Nx_core + 2; 
Ny = Ny_core + 2; 

% Reconstruct core and peripheral cell dimensions.
W_arr = ones(1, Nx) * W_core; W_arr(1) = W_edge; W_arr(end) = W_edge;
H_arr = ones(1, Ny) * H_core; H_arr(1) = H_edge; H_arr(end) = H_edge;

Base_P = 0.001; % Reference source density: 1 mW/um^2
Num_Pixels = Nx * Ny;
Num_Groups = 6; % Six virtual heat-source groups

% Times at which the response is sampled.
Target_Times = 1E-6:1E-6:13E-6;
Num_Times = length(Target_Times);

% --- 2. Load the solved response model ---
% The model name must match the output of Step1a.
ModelFile = sprintf('ICM_Hetero_6G_%dx%d.mph', Nx, Ny);
OutputFile = sprintf('G_Matrix_Hetero_%dx%d.mat', Nx, Ny);

fprintf('Load the response model: %s ...\n', ModelFile);
if ~exist(ModelFile, 'file')
    error('Model file %s is missing; run Step1a first.', ModelFile);
end
model = mphload(ModelFile);
fprintf('Response model loaded.\n');

% --- 3. Center and transverse-boundary sampling coordinates ---
X_bl = zeros(1, Nx); X_bl(1) = -W_edge;
for i = 2:Nx, X_bl(i) = X_bl(i-1) + W_arr(i-1); end

Y_bl = zeros(1, Ny); Y_bl(1) = -H_edge;
for j = 2:Ny, Y_bl(j) = Y_bl(j-1) + H_arr(j-1); end

Coords = zeros(3, Num_Pixels * 3);
idx = 0;
for j = 1:Ny
    for i = 1:Nx
        xc = X_bl(i) + W_arr(i)/2;
        yc = Y_bl(j) + H_arr(j)/2;
        % Sample the center and two x-boundaries at the silicon mid-plane (z = 110 nm).
        offsets_x = [0, -W_arr(i)/2, W_arr(i)/2]; 
        for k = 1:3
            idx = idx + 1;
            Coords(:, idx) = [xc + offsets_x(k); yc; 110];
        end
    end
end

% --- 4. Vectorized temperature extraction ---
fprintf('Extract %d temperature samples from dset2.\n', size(Coords,2));
tic;
% Select all outer solutions of the source-group parameter sweep.
Raw_Data = mphinterp(model, 'T', 'coord', Coords, ...
    'dataset', 'dset2', 't', Target_Times, 'outersolnum', 'all');
fprintf('Extraction time: %.2f s.\n', toc);

% --- 5. Reshape the response data ---
% The final dimension contains the six source-group responses.
G_4D = reshape((Raw_Data - T_amb) / Base_P, Num_Times, 3, Num_Pixels, Num_Groups);
G_Final = permute(G_4D, [2, 3, 1, 4]);

% --- 6. Package sampled responses ---
G_Data_Series = struct();
masks = generate_group_masks(Nx, Ny);
% Optimization points lie in the core region; exclude the peripheral cells.
monitor_indices = find(masks{1}); 

for k = 1:Num_Times
    G_Data_Series(k).time = Target_Times(k);
    
    % Sample 1: cell center.
    G_c_full = squeeze(G_Final(1, :, k, :)); 
    % Samples 2 and 3: the mean of the two transverse boundaries.
    G_b_full = squeeze(mean(G_Final(2:3, :, k, :), 1));
    
    % Responses at target-region cells for optimization.
    G_Data_Series(k).G_center = G_c_full(monitor_indices, :);
    G_Data_Series(k).G_base   = G_b_full(monitor_indices, :);
    
    % Responses at all cells for reconstruction and analysis.
    G_Data_Series(k).G_center_full = G_c_full;
    G_Data_Series(k).G_base_full   = G_b_full;
end

% --- 7. Save ---
save(OutputFile, 'G_Data_Series', 'masks', 'monitor_indices', 'T_amb', 'W_arr', 'H_arr');
fprintf('Response matrices saved: %s\n', OutputFile);
end

% --- Six-group masks ---
function masks = generate_group_masks(Nx, Ny)
    masks = cell(1, 6);
    
    % 1. Center
    mask_c = false(Nx, Ny); 
    mask_c(2:Nx-1, 2:Ny-1) = true; 
    masks{1} = mask_c;
    
    % 2. Top edge
    mask_top = false(Nx, Ny); 
    mask_top(2:Nx-1, Ny) = true; 
    masks{2} = mask_top;
    
    % 3. Bottom edge
    mask_bottom = false(Nx, Ny); 
    mask_bottom(2:Nx-1, 1) = true; 
    masks{3} = mask_bottom;
    
    % 4. Left and right edges
    mask_lr = false(Nx, Ny); 
    mask_lr([1, Nx], 2:Ny-1) = true; 
    masks{4} = mask_lr;
    
    % 5. Top corners
    mask_top_corner = false(Nx, Ny); 
    mask_top_corner([1, Nx], Ny) = true; 
    masks{5} = mask_top_corner;
    
    % 6. Bottom corners
    mask_bot_corner = false(Nx, Ny); 
    mask_bot_corner([1, Nx], 1) = true; 
    masks{6} = mask_bot_corner;
end

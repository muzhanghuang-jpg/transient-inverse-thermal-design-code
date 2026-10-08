function Step1b_Extract_G_Matrices_6G_Decoupled()
%STEP1B_EXTRACT_G_MATRICES_6G_DECOUPLED Build G matrices along the gradient axis.
%
% Fast mode expects Step1a to solve 6 cases at the gradient-axis extrema and
% centre:
%   Y mode: bottom / center / top Y-row groups (3 blocks × interior+side).
%   X mode: left / center / right X-column groups (3 blocks × interior+side).
% Edge blocks (first / last primary block) stay direct. Interior blocks are
% reconstructed from the centre block via translation along the gradient axis.
% Extended interpolation along that axis avoids truncating the translated
% response near boundaries.

import com.comsol.model.*
import com.comsol.model.util.*

powerWidthSettings = get_power_width_settings();
T_amb = powerWidthSettings.Ambient_Temp_K;
Nx_core = 31;
Ny_core = 21;
W_core = 2000;
H_core = 2000;
W_edge = 4000;
H_edge = 4000;
Translation_Fill_Mode = 'zero';

Nx = Nx_core + 2;
Ny = Ny_core + 2;

targetSettings = get_thermal_target_settings();
GradientDirection = targetSettings.GradientDirection;
isXMode = strcmpi(GradientDirection, 'X');

if isXMode
    Ndim_primary = Nx;
else
    Ndim_primary = Ny;
end
Center_Index = ceil(Ndim_primary / 2);
Center_Group_Indices = [2 * Center_Index - 1, 2 * Center_Index];
Edge_And_Center_Group_Indices = [1, 2, ...
    Center_Group_Indices, 2 * Ndim_primary - 1, 2 * Ndim_primary];

W_arr = ones(1, Nx) * W_core;
W_arr(1) = W_edge;
W_arr(end) = W_edge;

H_arr = ones(1, Ny) * H_core;
H_arr(1) = H_edge;
H_arr(end) = H_edge;

Base_P = 0.001; % W/um^2 probe density
Num_Pixels = Nx * Ny;
[masks, group_info, monitor_mask] = build_yrow_group_masks(Nx, Ny, GradientDirection);
Num_Groups = numel(masks);
monitor_indices = find(monitor_mask);

Target_Times = targetSettings.Time_Step_s:targetSettings.Time_Step_s:targetSettings.Target_Time_s;
Num_Times = numel(Target_Times);

ModelFile = sprintf('ICM_Hetero_YRow_%dx%d.mph', Nx, Ny);
OutputFile = sprintf('G_Matrix_Hetero_YRow_%dx%d.mat', Nx, Ny);

fprintf('Loading base model: %s (gradient direction: %s)\n', ModelFile, GradientDirection);
if ~exist(ModelFile, 'file')
    error('Cannot find %s. Run Step1a first.', ModelFile);
end
model = mphload(ModelFile);
fprintf('Base model loaded. Full matrix has %d groups along %s axis.\n', ...
    Num_Groups, GradientDirection);
fprintf('Fast edge+center probe groups are [%s].\n', ...
    num2str(Edge_And_Center_Group_Indices));

[X_centers, Y_centers] = build_pixel_centers(W_arr, H_arr);
CoordsActual = build_probe_coords(X_centers, Y_centers, W_arr);

fprintf('Interpolating actual array coordinates: %d points, %d time samples.\n', ...
    size(CoordsActual, 2), Num_Times);
tic;
Raw_Actual = mphinterp(model, 'T', 'coord', CoordsActual, ...
    'dataset', 'dset2', 't', Target_Times, 'outersolnum', 'all');
fprintf('Actual-coordinate interpolation elapsed time: %.2f s.\n', toc);

[G_Actual, Num_Source_Groups] = reshape_probe_data(Raw_Actual, T_amb, Base_P, ...
    Num_Times, 3, Num_Pixels);

if Num_Source_Groups == Num_Groups
    G_Matrix_Build_Mode = 'direct_full_sweep';
    G_Source_Group_Indices = 1:Num_Groups;
    fprintf('Detected full %d-group parametric sweep; no translation needed.\n', ...
        Num_Groups);
elseif Num_Source_Groups == numel(Edge_And_Center_Group_Indices)
    if isXMode
        G_Matrix_Build_Mode = 'translated_center_col_extended_x_with_direct_edges';
        fprintf(['Detected %d solved left/center/right-column groups; left and right ' ...
            'columns stay direct, interior columns use center-column translation.\n'], ...
            Num_Source_Groups);
    else
        G_Matrix_Build_Mode = 'translated_center_row_extended_y_with_direct_edges';
        fprintf(['Detected %d solved bottom/center/top-row groups; bottom and top ' ...
            'rows stay direct, interior rows use center-row translation.\n'], ...
            Num_Source_Groups);
    end
    G_Source_Group_Indices = Edge_And_Center_Group_Indices;
elseif Num_Source_Groups == numel(Center_Group_Indices)
    if isXMode
        G_Matrix_Build_Mode = 'translated_center_col_extended_x';
    else
        G_Matrix_Build_Mode = 'translated_center_row_extended_y';
    end
    G_Source_Group_Indices = Center_Group_Indices;
    fprintf(['Detected %d solved center-block groups; all blocks are reconstructed ' ...
        'from center-block translation.\n'], Num_Source_Groups);
else
    error(['Expected %d solved groups, edge+center groups [%s], or ' ...
        'center groups [%s], but mphinterp returned %d solved group(s).'], ...
        Num_Groups, num2str(Edge_And_Center_Group_Indices), ...
        num2str(Center_Group_Indices), Num_Source_Groups);
end

if strcmp(G_Matrix_Build_Mode, 'direct_full_sweep')
    G_Extended = [];
    Extended_Centers = [];
    Num_Extended_Probe_Pixels = 0;
else
    if isXMode
        Extended_Centers = build_extended_centers(X_centers, Center_Index);
        CoordsExtended = build_probe_coords_extended_x(Extended_Centers, Y_centers, W_core);
        Num_Extended_Probe_Pixels = numel(Extended_Centers) * Ny;
        axisLabel = 'col';
    else
        Extended_Centers = build_extended_centers(Y_centers, Center_Index);
        CoordsExtended = build_probe_coords(X_centers, Extended_Centers, W_arr);
        Num_Extended_Probe_Pixels = Nx * numel(Extended_Centers);
        axisLabel = 'row';
    end

    fprintf(['Interpolating extended center-%s coordinate range: %d samples ' ...
        'from %.0f nm to %.0f nm.\n'], axisLabel, ...
        numel(Extended_Centers), min(Extended_Centers), max(Extended_Centers));
    tic;
    Raw_Extended = mphinterp(model, 'T', 'coord', CoordsExtended, ...
        'dataset', 'dset2', 't', Target_Times, 'outersolnum', 'all');
    fprintf('Extended-coordinate interpolation elapsed time: %.2f s.\n', toc);

    [G_Extended, Num_Extended_Source_Groups] = reshape_probe_data( ...
        Raw_Extended, T_amb, Base_P, Num_Times, 3, Num_Extended_Probe_Pixels);

    if Num_Extended_Source_Groups ~= Num_Source_Groups
        error(['Actual interpolation returned %d solved groups, but extended ' ...
            'interpolation returned %d.'], Num_Source_Groups, ...
            Num_Extended_Source_Groups);
    end
end

G_Data_Series = struct();
for tIdx = 1:Num_Times
    G_Data_Series(tIdx).time = Target_Times(tIdx);

    G_c_actual = reshape(G_Actual(1, :, tIdx, :), ...
        Num_Pixels, Num_Source_Groups);
    G_b_actual = reshape(mean(G_Actual(2:3, :, tIdx, :), 1), ...
        Num_Pixels, Num_Source_Groups);

    if strcmp(G_Matrix_Build_Mode, 'direct_full_sweep')
        G_c_full = G_c_actual;
        G_b_full = G_b_actual;
    else
        G_c_extended = reshape(G_Extended(1, :, tIdx, :), ...
            Num_Extended_Probe_Pixels, Num_Source_Groups);
        G_b_extended = reshape(mean(G_Extended(2:3, :, tIdx, :), 1), ...
            Num_Extended_Probe_Pixels, Num_Source_Groups);
        if isXMode
            G_c_full = expand_center_translated_matrix_x(G_c_actual, ...
                G_c_extended, G_Source_Group_Indices, masks, group_info, ...
                W_arr, H_arr, X_centers, Extended_Centers, Center_Index, ...
                Translation_Fill_Mode);
            G_b_full = expand_center_translated_matrix_x(G_b_actual, ...
                G_b_extended, G_Source_Group_Indices, masks, group_info, ...
                W_arr, H_arr, X_centers, Extended_Centers, Center_Index, ...
                Translation_Fill_Mode);
        else
            G_c_full = expand_center_translated_matrix_y(G_c_actual, ...
                G_c_extended, G_Source_Group_Indices, masks, group_info, ...
                W_arr, H_arr, Y_centers, Extended_Centers, Center_Index, ...
                Translation_Fill_Mode);
            G_b_full = expand_center_translated_matrix_y(G_b_actual, ...
                G_b_extended, G_Source_Group_Indices, masks, group_info, ...
                W_arr, H_arr, Y_centers, Extended_Centers, Center_Index, ...
                Translation_Fill_Mode);
        end
    end

    G_Data_Series(tIdx).G_center = G_c_full(monitor_indices, :);
    G_Data_Series(tIdx).G_base = G_b_full(monitor_indices, :);
    G_Data_Series(tIdx).G_center_full = G_c_full;
    G_Data_Series(tIdx).G_base_full = G_b_full;
end

save(OutputFile, 'G_Data_Series', 'masks', 'group_info', ...
    'monitor_indices', 'T_amb', 'W_arr', 'H_arr', ...
    'G_Matrix_Build_Mode', 'G_Source_Group_Indices', 'Center_Index', ...
    'Translation_Fill_Mode', 'Extended_Centers', 'GradientDirection');
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
    % Loop order (j outer, i inner) preserves the (Nx × Ny) reshape used by
    % expand_center_translated_matrix_y. Each pixel column i uses its own
    % W_arr(i) for the left/right X offsets so edge columns are sampled
    % correctly.
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

function Coords = build_probe_coords_extended_x(Extended_X_Centers, Y_centers, W_default)
    % X mode extended probe. Loop order (xIdx outer, j inner) so reshape to
    % (Ny × NumExtX) gives sourceMapExtended(j, xIdx). At extended X positions
    % no actual pixel column exists, so a constant W_default (= W_core, the
    % center column's width) is used for the left/right X offsets. Edge X
    % columns are probed directly via the actual-coordinate path, where the
    % per-column W_arr(i) is honoured.
    NxExt = numel(Extended_X_Centers);
    Ny = numel(Y_centers);
    Coords = zeros(3, NxExt * Ny * 3);

    idx = 0;
    for xIdx = 1:NxExt
        for j = 1:Ny
            offsets_x = [0, -W_default / 2, W_default / 2];
            for k = 1:3
                idx = idx + 1;
                Coords(:, idx) = [Extended_X_Centers(xIdx) + offsets_x(k); ...
                    Y_centers(j); 110];
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

function Extended_Centers = build_extended_centers(centers, centerIdx)
    % Extended primary-axis coordinates required by the translation lookup:
    %   queryC = centers(obs) - centers(target) + centers(centerIdx)
    % over all (obs, target) pairs.
    N = numel(centers);
    values = [];

    for targetBlock = 1:N
        for obsBlock = 1:N
            values(end + 1) = centers(obsBlock) - ...
                centers(targetBlock) + centers(centerIdx); %#ok<AGROW>
        end
    end

    Extended_Centers = unique(values);
end

function G_full = expand_center_translated_matrix_y(G_actual, G_extended, ...
    source_group_indices, masks, group_info, W_arr, H_arr, Y_centers, ...
    Extended_Y_Centers, centerIdx, fillMode)

    Nx = numel(W_arr);
    Ny = numel(H_arr);
    Num_Pixels = Nx * Ny;
    Num_Groups = numel(masks);

    areaMatrix = (W_arr(:) * H_arr(:)') / (1000 ^ 2);
    centerInteriorGroup = 2 * centerIdx - 1;
    centerSideGroup = 2 * centerIdx;
    centerInteriorCol = find(source_group_indices == centerInteriorGroup, 1);
    centerSideCol = find(source_group_indices == centerSideGroup, 1);

    if isempty(centerInteriorCol) || isempty(centerSideCol)
        error('Center-row source columns [%d %d] are required for Y translation.', ...
            centerInteriorGroup, centerSideGroup);
    end

    G_full = zeros(Num_Pixels, Num_Groups);
    for groupIdx = 1:Num_Groups
        directSourceCol = find(source_group_indices == groupIdx, 1);
        if ~isempty(directSourceCol)
            G_full(:, groupIdx) = G_actual(:, directSourceCol);
            continue;
        end

        groupType = group_info(groupIdx).type;
        targetRow = group_info(groupIdx).row_j;
        targetArea = sum(areaMatrix(masks{groupIdx}));

        if strcmp(groupType, 'interior')
            sourceCol = centerInteriorCol;
            sourceArea = sum(areaMatrix(masks{centerInteriorGroup}));
        elseif strcmp(groupType, 'side_lr')
            sourceCol = centerSideCol;
            sourceArea = sum(areaMatrix(masks{centerSideGroup}));
        else
            error('Unsupported Y-row group type: %s', groupType);
        end

        sourceMapExtended = reshape(G_extended(:, sourceCol), ...
            Nx, numel(Extended_Y_Centers));
        translatedMap = zeros(Nx, Ny);

        for obsRow = 1:Ny
            queryY = Y_centers(obsRow) - Y_centers(targetRow) + ...
                Y_centers(centerIdx);
            yIdx = find(Extended_Y_Centers == queryY, 1);

            if ~isempty(yIdx)
                translatedMap(:, obsRow) = sourceMapExtended(:, yIdx) * ...
                    (targetArea / sourceArea);
            elseif strcmp(fillMode, 'nearest')
                [~, yIdx] = min(abs(Extended_Y_Centers - queryY));
                translatedMap(:, obsRow) = sourceMapExtended(:, yIdx) * ...
                    (targetArea / sourceArea);
            elseif ~strcmp(fillMode, 'zero')
                error('Unsupported translation fill mode: %s', fillMode);
            end
        end

        G_full(:, groupIdx) = translatedMap(:);
    end
end

function G_full = expand_center_translated_matrix_x(G_actual, G_extended, ...
    source_group_indices, masks, group_info, W_arr, H_arr, X_centers, ...
    Extended_X_Centers, centerIdx, fillMode)
    % X mode translation: centre column's response is shifted along X to fill
    % every interior column's G column. sourceMapExtended is (Ny × NumExtX):
    % rows are pixel Y indices, columns are extended X positions.

    Nx = numel(W_arr);
    Ny = numel(H_arr);
    Num_Pixels = Nx * Ny;
    Num_Groups = numel(masks);

    areaMatrix = (W_arr(:) * H_arr(:)') / (1000 ^ 2);
    centerInteriorGroup = 2 * centerIdx - 1;
    centerSideGroup = 2 * centerIdx;
    centerInteriorCol = find(source_group_indices == centerInteriorGroup, 1);
    centerSideCol = find(source_group_indices == centerSideGroup, 1);

    if isempty(centerInteriorCol) || isempty(centerSideCol)
        error('Center-column source columns [%d %d] are required for X translation.', ...
            centerInteriorGroup, centerSideGroup);
    end

    G_full = zeros(Num_Pixels, Num_Groups);
    for groupIdx = 1:Num_Groups
        directSourceCol = find(source_group_indices == groupIdx, 1);
        if ~isempty(directSourceCol)
            G_full(:, groupIdx) = G_actual(:, directSourceCol);
            continue;
        end

        groupType = group_info(groupIdx).type;
        targetCol = group_info(groupIdx).col_i;
        targetArea = sum(areaMatrix(masks{groupIdx}));

        if strcmp(groupType, 'interior')
            sourceCol = centerInteriorCol;
            sourceArea = sum(areaMatrix(masks{centerInteriorGroup}));
        elseif strcmp(groupType, 'side_tb')
            sourceCol = centerSideCol;
            sourceArea = sum(areaMatrix(masks{centerSideGroup}));
        else
            error('Unsupported X-column group type: %s', groupType);
        end

        % Probe order in build_probe_coords_extended_x: (xIdx outer, j inner).
        % MATLAB column-major linearisation puts j as the fastest-varying
        % index, so reshape to (Ny × NumExtX) reads back as (j, xIdx).
        sourceMapExtended = reshape(G_extended(:, sourceCol), ...
            Ny, numel(Extended_X_Centers));
        translatedMap = zeros(Nx, Ny);

        for obsCol = 1:Nx
            queryX = X_centers(obsCol) - X_centers(targetCol) + ...
                X_centers(centerIdx);
            xIdx = find(Extended_X_Centers == queryX, 1);

            if ~isempty(xIdx)
                translatedMap(obsCol, :) = sourceMapExtended(:, xIdx).' * ...
                    (targetArea / sourceArea);
            elseif strcmp(fillMode, 'nearest')
                [~, xIdx] = min(abs(Extended_X_Centers - queryX));
                translatedMap(obsCol, :) = sourceMapExtended(:, xIdx).' * ...
                    (targetArea / sourceArea);
            elseif ~strcmp(fillMode, 'zero')
                error('Unsupported translation fill mode: %s', fillMode);
            end
        end

        G_full(:, groupIdx) = translatedMap(:);
    end
end

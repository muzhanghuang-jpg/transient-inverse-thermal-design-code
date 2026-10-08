function [masks, group_info, monitor_mask] = build_yrow_group_masks(Nx, Ny, gradientDirection)
%BUILD_YROW_GROUP_MASKS Control-group masks aligned with the gradient axis.
%
% gradientDirection: 'X' or 'Y' (defaults to 'Y' for backward compatibility).
%   'Y' mode: 2*Ny groups, each Y row split into interior + LR-edge.
%   'X' mode: 2*Nx groups, each X column split into interior + TB-edge.
%
% Group ordering within either mode:
%   2*k-1 -> interior of primary block k (avoids edge in the secondary axis)
%   2*k   -> edge of primary block k (the two pixels at the secondary-axis ends)
%
% group_info fields:
%   index        - group number (1..numGroups)
%   primary_idx  - row j (Y mode) or column i (X mode), the primary index
%   type         - 'interior' or 'side_lr' (Y mode) / 'side_tb' (X mode)
%   row_j, col_i - filled only on the matching axis for human-readable debugging

    if nargin < 3 || isempty(gradientDirection)
        gradientDirection = 'Y';
    end

    isXMode = strcmpi(gradientDirection, 'X');
    if isXMode
        Ndim = Nx;
        edgeTypeStr = 'side_tb';
    else
        Ndim = Ny;
        edgeTypeStr = 'side_lr';
    end

    numGroups = 2 * Ndim;
    masks = cell(1, numGroups);
    group_info = repmat(struct('index', 0, 'primary_idx', 0, 'type', '', ...
        'row_j', 0, 'col_i', 0), 1, numGroups);

    for k_dim = 1:Ndim
        interiorIndex = 2 * k_dim - 1;
        sideIndex = 2 * k_dim;
        interiorMask = false(Nx, Ny);
        sideMask = false(Nx, Ny);

        if isXMode
            interiorMask(k_dim, 2:Ny-1) = true;
            sideMask(k_dim, [1, Ny]) = true;
            group_info(interiorIndex).col_i = k_dim;
            group_info(sideIndex).col_i = k_dim;
        else
            interiorMask(2:Nx-1, k_dim) = true;
            sideMask([1, Nx], k_dim) = true;
            group_info(interiorIndex).row_j = k_dim;
            group_info(sideIndex).row_j = k_dim;
        end

        masks{interiorIndex} = interiorMask;
        masks{sideIndex} = sideMask;

        group_info(interiorIndex).index = interiorIndex;
        group_info(interiorIndex).type = 'interior';
        group_info(interiorIndex).primary_idx = k_dim;

        group_info(sideIndex).index = sideIndex;
        group_info(sideIndex).type = edgeTypeStr;
        group_info(sideIndex).primary_idx = k_dim;
    end

    monitor_mask = false(Nx, Ny);
    monitor_mask(2:Nx-1, 2:Ny-1) = true;
end

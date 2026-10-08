function [masks, group_info, monitor_mask] = build_yrow_group_masks(Nx, Ny)
%BUILD_YROW_GROUP_MASKS Symmetric 4-group control map for the heterogeneous array.
%
% Group order:
%   1: edge_interior  - top+bottom row interior pixels merged (i=2..Nx-1, j in {1,Ny})
%   2: edge_side      - top+bottom row L/R-edge pixels merged (i in {1,Nx}, j in {1,Ny})
%   3: middle_interior- middle ROI block (i=2..Nx-1, j=2..Ny-1)
%   4: middle_side    - middle L/R-edge column (i in {1,Nx}, j=2..Ny-1)

    numGroups = 4;
    masks = cell(1, numGroups);
    group_info = repmat(struct('index', 0, 'type', ''), 1, numGroups);

    edgeInterior = false(Nx, Ny);
    edgeInterior(2:Nx-1, [1, Ny]) = true;
    masks{1} = edgeInterior;
    group_info(1).index = 1;
    group_info(1).type = 'edge_interior';

    edgeSide = false(Nx, Ny);
    edgeSide([1, Nx], [1, Ny]) = true;
    masks{2} = edgeSide;
    group_info(2).index = 2;
    group_info(2).type = 'edge_side';

    middleInterior = false(Nx, Ny);
    middleInterior(2:Nx-1, 2:Ny-1) = true;
    masks{3} = middleInterior;
    group_info(3).index = 3;
    group_info(3).type = 'middle_interior';

    middleSide = false(Nx, Ny);
    middleSide([1, Nx], 2:Ny-1) = true;
    masks{4} = middleSide;
    group_info(4).index = 4;
    group_info(4).type = 'middle_side';

    monitor_mask = false(Nx, Ny);
    monitor_mask(2:Nx-1, 2:Ny-1) = true;
end

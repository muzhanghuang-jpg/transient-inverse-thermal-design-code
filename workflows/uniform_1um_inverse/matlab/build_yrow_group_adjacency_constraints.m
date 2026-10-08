function Constraints_List = build_yrow_group_adjacency_constraints()
%BUILD_YROW_GROUP_ADJACENCY_CONSTRAINTS Contrast pairs for the 4-group control map.
%
% Constraints are returned as directed pairs [a b], meaning P(a) <= C * P(b).
% Both directions are included for each local adjacency.
%
% Group layout (see build_yrow_group_masks):
%   1: edge_interior, 2: edge_side, 3: middle_interior, 4: middle_side

    pairs = [1 2; 3 4; 1 3; 2 4];

    Constraints_List = [];
    for k = 1:size(pairs, 1)
        a = pairs(k, 1);
        b = pairs(k, 2);
        Constraints_List = [Constraints_List; a, b; b, a]; %#ok<AGROW>
    end

    Constraints_List = unique(Constraints_List, 'rows');
end

function Constraints_List = build_yrow_group_adjacency_constraints(numBlocks)
%BUILD_YROW_GROUP_ADJACENCY_CONSTRAINTS Contrast constraints for control groups.
%
% numBlocks: number of primary-direction blocks. Pass Ny in Y-gradient mode,
%   Nx in X-gradient mode. The adjacency structure is direction-agnostic:
%   interior <-> side within each block, plus interior <-> interior and
%   side <-> side between adjacent blocks.
%
% Constraints are returned as directed pairs [a b], meaning P(a) <= C * P(b).
% Both directions are included for each local adjacency.

    Constraints_List = [];

    for k = 1:numBlocks
        interior = 2 * k - 1;
        side = 2 * k;
        Constraints_List = add_directed_pair(Constraints_List, interior, side);

        if k < numBlocks
            nextInterior = 2 * (k + 1) - 1;
            nextSide = 2 * (k + 1);
            Constraints_List = add_directed_pair(Constraints_List, interior, nextInterior);
            Constraints_List = add_directed_pair(Constraints_List, side, nextSide);
        end
    end

    Constraints_List = unique(Constraints_List, 'rows');
end

function pairs = add_directed_pair(pairs, a, b)
    pairs = [pairs; a, b; b, a];
end

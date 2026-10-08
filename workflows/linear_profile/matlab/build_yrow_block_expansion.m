function [E, block_info, row_to_block] = build_yrow_block_expansion(Ny, Y_Row_Block_Size)
%BUILD_YROW_BLOCK_EXPANSION Map block-level Y controls to per-row groups.
%
% Full group order is inherited from build_yrow_group_masks:
%   2*j - 1: interior group on Y row j
%   2*j    : merged left/right side group on Y row j
%
% Block group order:
%   2*b - 1: interior group for Y-row block b
%   2*b    : merged left/right side group for Y-row block b

    if nargin < 2 || isempty(Y_Row_Block_Size)
        Y_Row_Block_Size = 1;
    end
    if ~(isscalar(Y_Row_Block_Size) && isfinite(Y_Row_Block_Size) && ...
            Y_Row_Block_Size >= 1 && Y_Row_Block_Size == round(Y_Row_Block_Size))
        error('Y_Row_Block_Size must be a positive integer.');
    end

    numFullGroups = 2 * Ny;
    numBlocks = ceil(Ny / Y_Row_Block_Size);
    numBlockGroups = 2 * numBlocks;

    E = zeros(numFullGroups, numBlockGroups);
    row_to_block = zeros(1, Ny);
    block_info = repmat(struct( ...
        'block_index', 0, ...
        'start_row_j', 0, ...
        'end_row_j', 0, ...
        'interior_group', 0, ...
        'side_lr_group', 0), 1, numBlocks);

    for b = 1:numBlocks
        startRow = (b - 1) * Y_Row_Block_Size + 1;
        endRow = min(b * Y_Row_Block_Size, Ny);
        interiorBlockGroup = 2 * b - 1;
        sideBlockGroup = 2 * b;

        block_info(b).block_index = b;
        block_info(b).start_row_j = startRow;
        block_info(b).end_row_j = endRow;
        block_info(b).interior_group = interiorBlockGroup;
        block_info(b).side_lr_group = sideBlockGroup;

        for j = startRow:endRow
            row_to_block(j) = b;
            E(2 * j - 1, interiorBlockGroup) = 1;
            E(2 * j, sideBlockGroup) = 1;
        end
    end
end

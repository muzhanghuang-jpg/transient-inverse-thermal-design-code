function D = build_column_smoothness_matrix(Nx)
%BUILD_COLUMN_SMOOTHNESS_MATRIX Second-difference regulariser for X-mode
% interior groups AND TB-edge groups across interior columns
% (col 2..Nx-1). Edge columns col 1 and col Nx are excluded entirely.
%
% Returns a (2*(Nx-4)) x (2*Nx) matrix D. Multiplied by the X-mode group
% power vector P (group ordering: 2*i-1 = column i interior, 2*i =
% column i TB-edge), the rows evaluate the discrete Laplacian along X:
%
%   Rows 1..Nx-4:               interior 2nd-difference
%       D(k, 2*(k+1)-1)   = +1     interior of col k+1
%       D(k, 2*(k+2)-1)   = -2     interior of col k+2 (centre)
%       D(k, 2*(k+3)-1)   = +1     interior of col k+3
%
%   Rows Nx-3..2*(Nx-4):        side (TB-edge) 2nd-difference
%       D(Nx-4+k, 2*(k+1)) = +1    TB-edge of col k+1
%       D(Nx-4+k, 2*(k+2)) = -2    TB-edge of col k+2 (centre)
%       D(Nx-4+k, 2*(k+3)) = +1    TB-edge of col k+3
%
% Both blocks penalise |P(i-1) - 2*P(i) + P(i+1)|, allowing P to vary
% linearly with X (consistent with the X-gradient ramp) while suppressing
% the alternating-column zigzag and higher-order curvature noise.
%
% Excluded from D (unconstrained by the regulariser):
%   - Edge columns col 1 and col Nx (groups 1, 2, 2*Nx-1, 2*Nx). They
%     compensate X-boundary heat loss and are expected to carry power
%     densities sharply different from their interior neighbours. The
%     side groups of col 1 and col Nx are also the corner pixels
%     (W_edge x H_edge), structurally different from the TB-edge pixels
%     of interior columns (W_core x H_edge).
%
% Use only in X mode. Y mode does not call this helper.

    if ~(isscalar(Nx) && isfinite(Nx) && Nx >= 5 && Nx == floor(Nx))
        error('Nx must be an integer >= 5 to admit a second-difference row.');
    end

    Num_Groups = 2 * Nx;
    rows_per_block = Nx - 4;
    num_rows = 2 * rows_per_block;
    D = zeros(num_rows, Num_Groups);

    for k = 1:rows_per_block
        % Interior second-difference (top half of D)
        D(k, 2*(k+1) - 1) = +1;
        D(k, 2*(k+2) - 1) = -2;
        D(k, 2*(k+3) - 1) = +1;

        % Side (TB-edge) second-difference (bottom half of D)
        D(rows_per_block + k, 2*(k+1)) = +1;
        D(rows_per_block + k, 2*(k+2)) = -2;
        D(rows_per_block + k, 2*(k+3)) = +1;
    end
end

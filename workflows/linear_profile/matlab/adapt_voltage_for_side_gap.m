function [Voltage_Use, width_map_limited, adjusted] = adapt_voltage_for_side_gap( ...
    width_map_raw, W_arr, Voltage_Use, Min_Wire_Width_nm, Min_Side_Gap_nm, context)
%ADAPT_VOLTAGE_FOR_SIDE_GAP Increase voltage if raw widths would erase side gaps.

if nargin < 6
    context = 'wire width generation';
end

if Voltage_Use <= 0
    error('Voltage_Use must be positive to adapt wire widths.');
end

W_max_by_row = max(Min_Wire_Width_nm, W_arr(:) - 2 * Min_Side_Gap_nm);
W_max_map = repmat(W_max_by_row, 1, size(width_map_raw, 2));

positive_width = width_map_raw > 0;
over_max = positive_width & (width_map_raw > W_max_map);
adjusted = any(over_max(:));

if adjusted
    % For fixed target power, wire width scales as 1/V^2, so sqrt(max ratio)
    % raises voltage just enough to bring the widest cell down to its limit.
    voltage_scale = sqrt(max(width_map_raw(over_max) ./ W_max_map(over_max)));
    old_voltage = Voltage_Use;
    Voltage_Use = Voltage_Use * voltage_scale;
    width_map_raw(positive_width) = width_map_raw(positive_width) / (voltage_scale ^ 2);

    fprintf(['%s: increased Voltage_Use from %.3f V to %.3f V to preserve ' ...
        '%.0f nm side gaps (%d cell(s) exceeded max width).\n'], ...
        context, old_voltage, Voltage_Use, Min_Side_Gap_nm, nnz(over_max));
end

width_map_limited = min(max(width_map_raw, Min_Wire_Width_nm), W_max_map);

end

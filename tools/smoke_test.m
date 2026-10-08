function smoke_test()
root = fileparts(fileparts(mfilename('fullpath')));
files = dir(fullfile(root, '**', '*.m'));
errors = {};
for k = 1:numel(files)
    file = fullfile(files(k).folder,files(k).name);
    findings = checkcode(file, '-id');
    for j = 1:numel(findings)
        if startsWith(findings(j).id, 'SYNER') || strcmp(findings(j).id, 'ENDCT')
            errors{end+1} = sprintf('%s:%d %s', file, findings(j).line, findings(j).message); %#ok<AGROW>
        end
    end
end
assert(isempty(errors), strjoin(errors,newline));
previous = path;
cleanup = onCleanup(@() path(previous));
square = fullfile(root,'workflows','uniform_1um_inverse','matlab');
addpath(square);
s = get_thermal_target_settings();
assert(s.Base_Temp_K == 900 && s.Target_Time_s == 13e-6 && s.Y_Gradient_K_per_um == 0);
rmpath(square);
clear get_thermal_target_settings
addpath(fullfile(root,'workflows','linear_profile','matlab'));
s = get_thermal_target_settings();
assert(strcmp(s.GradientDirection,'X') && s.Y_Gradient_K_per_um == 2 && s.Target_Time_s == 13e-6);
fprintf('PASS: %d MATLAB files parsed; square/gradient target settings checked. No FEM solve performed.\n',numel(files));
clear cleanup
end

function pathValue = comsol_mli_path()
% Resolve the installation without assuming the office computer's path.
pathValue = getenv('COMSOL_MLI_PATH');
if isempty(pathValue)
    located = which('mphstart');
    if ~isempty(located), pathValue = fileparts(located); end
end
if isempty(pathValue) && ispc
    candidate = fullfile(getenv('ProgramFiles'), 'COMSOL', 'COMSOL64', 'Multiphysics', 'mli');
    if isfolder(candidate), pathValue = candidate; end
end
if isempty(pathValue) || ~isfolder(pathValue)
    error('thermal_design:MissingLiveLink', 'Set COMSOL_MLI_PATH to your LiveLink mli folder.');
end
end

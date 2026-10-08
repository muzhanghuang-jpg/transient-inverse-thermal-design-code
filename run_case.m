function run_case(caseId)
% Run one isolated case with the original numerical settings.
root = fileparts(mfilename('fullpath'));
ids = {'uniform_1um_inverse','uniform_1um_plain','rectangular_2um_1d1r','linear_profile'};
caseId = validatestring(caseId, ids);
caseRoot = fullfile(root, 'workflows', caseId);
previousPath = path;
previousDir = pwd;
cleanup = onCleanup(@() restore_environment(previousPath, previousDir));
for k = 1:numel(ids)
    folder = fullfile(root, 'workflows', ids{k}, 'matlab');
    if contains([path pathsep], [folder pathsep]), rmpath(folder); end
end
addpath(fullfile(root, 'tools'), fullfile(caseRoot, 'matlab'));
if ismember(caseId, {'uniform_1um_inverse','linear_profile'})
    Run_Full_Auto_Workflow();
else
    addpath(comsol_mli_path());
    mphstart(comsol_server_port());
    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
    runDir = fullfile(caseRoot, 'results', ['Run_Full_' stamp]);
    assert(~isfolder(runDir), 'Output folder already exists.');
    mkdir(runDir);
    cd(runDir);
    diary(fullfile(runDir, 'workflow.log'));
    diaryCleanup = onCleanup(@() diary('off'));
    if strcmp(caseId, 'uniform_1um_plain')
        Build_Plain_Heater(20.92);
    else
        Step1a_Generate_ICM_Model_Hetero(130, 1, 1);
        Step1b_Extract_G_Matrices_6G_Decoupled();
        Step2b_optimization_Physical_WidthGen_Hetero();
        Step3_Build_Verify_Model_Hetero_diode(0, 130, 1, 1);
        Step4_Iterative_Correction_Hetero(0, 0.5);
        Run_Full_Auto_Workflow();
    end
    clear diaryCleanup
end
clear cleanup
end

function restore_environment(previousPath, previousDir)
cd(previousDir);
path(previousPath);
end

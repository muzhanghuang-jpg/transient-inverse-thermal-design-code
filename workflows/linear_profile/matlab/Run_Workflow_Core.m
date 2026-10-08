function Run_Workflow_Core(options)
%RUN_WORKFLOW_CORE Shared implementation for full and resumed workflows.

    if nargin < 1
        options = struct();
    end
    options = normalize_workflow_options(options);

    originalDir = pwd;
    scriptDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(scriptDir);

    addpath(scriptDir);
    ensure_comsol_livelink_connected();
    ensure_project_folders(projectRoot);

    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    [runDir, runName] = create_timestamped_run_dir(projectRoot, timestamp, options.RunPrefix);

    cleanupCd = onCleanup(@() cd(originalDir));
    cd(runDir);

    logFile = fullfile(runDir, [runName '.log']);
    diary(logFile);
    cleanupDiary = onCleanup(@() diary('off'));

    fprintf('\n==================================================\n');
    fprintf('COMSOL-MATLAB timestamped workflow started\n');
    fprintf('Mode: %s\n', options.Mode);
    fprintf('Project root: %s\n', projectRoot);
    fprintf('Run directory: %s\n', runDir);
    fprintf('Log file: %s\n', logFile);
    fprintf('==================================================\n');

    try
        config = create_workflow_config();

        if strcmp(options.Mode, 'resume')
            fprintf('\nResume source run directory: %s\n', options.SourceRunDir);
            fprintf('Resume start step: %s\n', options.StartStep);
            fprintf('Resume start iter: %d\n', options.StartIter);
            copy_resume_artifacts(options.SourceRunDir, runDir, options.StartStep, options.StartIter);
        end

        fprintf('\nMaterial parameters:\n');
        fprintf('  k_SOI_val    = %.6g W/(m*K)\n', config.k_SOI_val);
        fprintf('  k_Hdoped_val = %.6g\n', config.k_Hdoped_val);
        fprintf('  k_Ldoped_val = %.6g\n', config.k_Ldoped_val);
        fprintf('  Initial Gain = %.3f\n', config.Gain);
        fprintf('  All generated files for this run are written under: %s\n', runDir);

        startOrder = step_order(options.StartStep);

        if startOrder <= step_order('Step1a')
            fprintf('\n[1/4] Step1a: build and solve base heterogeneous thermal model\n');
            Step1a_Generate_ICM_Model_Hetero(config.k_SOI_val, ...
                config.k_Hdoped_val, config.k_Ldoped_val);
        else
            fprintf('\n[1/4] Step1a skipped; using copied base model artifacts.\n');
        end

        if startOrder <= step_order('Step1b')
            fprintf('\n[2/4] Step1b: extract control-group G matrix\n');
            Step1b_Extract_G_Matrices_6G_Decoupled();
        else
            fprintf('\n[2/4] Step1b skipped; using copied G matrix artifacts.\n');
        end

        if startOrder <= step_order('Step2b')
            fprintf('\n[3/4] Step2b: compute initial optimized power density and linewidth map\n');
            Step2b_optimization_Physical_WidthGen_Hetero();
        else
            fprintf('\n[3/4] Step2b skipped; using copied optimization/correction artifacts.\n');
        end

        fprintf('\n[4/4] Closed-loop verification starts from iter = %d\n', options.StartIter);
        run_closed_loop(config, options.StartIter, options.StartStep);

    catch ME
        fprintf(2, '\nWorkflow failed: %s\n', ME.message);
        for k = 1:numel(ME.stack)
            fprintf(2, '  at %s:%d\n', ME.stack(k).file, ME.stack(k).line);
        end
        rethrow(ME);
    end

    clear cleanupDiary cleanupCd
end

function options = normalize_workflow_options(options)
    if ~isfield(options, 'Mode') || isempty(options.Mode)
        options.Mode = 'full';
    end
    options.Mode = validatestring(options.Mode, {'full', 'resume'});

    if strcmp(options.Mode, 'full')
        options.StartStep = 'Step1a';
        options.StartIter = 0;
        options.SourceRunDir = '';
        options.RunPrefix = 'Run_Full_Auto_Workflow';
        return;
    end

    if ~isfield(options, 'SourceRunDir') || isempty(options.SourceRunDir)
        error('Resume mode requires SourceRunDir.');
    end
    sourceRunDir = char(options.SourceRunDir);
    if exist(sourceRunDir, 'dir') ~= 7
        error('Resume source run directory does not exist: %s', sourceRunDir);
    end
    options.SourceRunDir = char(java.io.File(sourceRunDir).getCanonicalPath());

    if ~isfield(options, 'StartStep') || isempty(options.StartStep)
        options.StartStep = 'Step2b';
    end
    options.StartStep = normalize_start_step(options.StartStep);

    if ~isfield(options, 'StartIter') || isempty(options.StartIter)
        options.StartIter = 0;
    end
    options.StartIter = validate_nonnegative_integer(options.StartIter, 'StartIter');
    options.RunPrefix = 'Run_Resume_Workflow';
end

function startStep = normalize_start_step(startStep)
    if isnumeric(startStep)
        order = validate_nonnegative_integer(startStep, 'startStep');
        names = {'Step1a', 'Step1b', 'Step2b', 'Step3', 'Step4'};
        if order < 1 || order > numel(names)
            error('Numeric startStep must be 1 through 5.');
        end
        startStep = names{order};
        return;
    end

    raw = lower(strtrim(char(startStep)));
    raw = regexprep(raw, '[\s_-]', '');
    switch raw
        case {'1', 'step1', 'step1a'}
            startStep = 'Step1a';
        case {'2', 'step1b'}
            startStep = 'Step1b';
        case {'3', 'step2', 'step2b'}
            startStep = 'Step2b';
        case {'4', 'step3'}
            startStep = 'Step3';
        case {'5', 'step4'}
            startStep = 'Step4';
        otherwise
            error('Unsupported startStep "%s". Use Step1a, Step1b, Step2b, Step3, or Step4.', startStep);
    end
end

function value = validate_nonnegative_integer(value, name)
    if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0 && value == floor(value))
        error('%s must be a nonnegative integer scalar.', name);
    end
end

function config = create_workflow_config()
    config.k_SOI_val = 130;
    config.k_Hdoped_val = 1;
    config.k_Ldoped_val = 1;

    config.Max_Iterations = 15;
    config.Target_RMSE_K = 3.0;
    config.Target_MaxAbsError_K = 5.0;
    config.Gain = 0.8;
    config.Min_Gain = 0.05;
    config.Max_Rollback_Attempts_Per_Iter = 6;
    config.Stall_Patience = 3;
    config.Min_Score_Improvement = 1e-3;
end

function run_closed_loop(config, startIter, startStep)
    iter = startIter;
    exitStatus = 'running';
    Gain = config.Gain;
    rollbackAttemptsThisIter = 0;
    stalledAcceptedIterations = 0;
    prevRMSE = inf;
    prevMaxAbsError = inf;

    lastAccepted.iter = -1;
    lastAccepted.rmse = inf;
    lastAccepted.maxAbsError = inf;
    lastAccepted.score = inf;
    lastAccepted.gain = Gain;
    bestAccepted = lastAccepted;

    firstLoopStartsAtStep4 = strcmp(startStep, 'Step4');

    while iter <= config.Max_Iterations
        fprintf('\n==================================================\n');
        fprintf('Evaluating iter = %d with Gain = %.3f\n', iter, Gain);
        fprintf('==================================================\n');

        if firstLoopStartsAtStep4
            fprintf('Step3 skipped for first resumed loop; using copied verification model for iter %d.\n', iter);
            firstLoopStartsAtStep4 = false;
        else
            Step3_Build_Verify_Model_Hetero_diode(iter, ...
                config.k_SOI_val, config.k_Hdoped_val, config.k_Ldoped_val);
        end

        [currentRMSE, maxPosError, maxNegError] = ...
            Step4_Iterative_Correction_Hetero(iter, Gain);

        currentMaxAbsError = max(abs([maxPosError, maxNegError]));

        fprintf('\nIter %d measured ROI error:\n', iter);
        fprintf('  RMSE              = %.4f K\n', currentRMSE);
        fprintf('  Max positive error = %+.4f K\n', maxPosError);
        fprintf('  Max negative error = %+.4f K\n', maxNegError);
        fprintf('  Max abs error      = %.4f K\n', currentMaxAbsError);

        if iter > startIter && currentMaxAbsError > prevMaxAbsError + 1e-9
            rollbackAttemptsThisIter = rollbackAttemptsThisIter + 1;
            fprintf('\nRollback triggered for iter %d:\n', iter);
            fprintf('  Current max abs error %.4f K is worse than previous accepted value %.4f K.\n', ...
                currentMaxAbsError, prevMaxAbsError);
            fprintf('  Rollback attempt for this iter: %d/%d.\n', ...
                rollbackAttemptsThisIter, config.Max_Rollback_Attempts_Per_Iter);

            if rollbackAttemptsThisIter >= config.Max_Rollback_Attempts_Per_Iter
                exitStatus = 'rollback_limit';
                fprintf(['  Stopping without accepting iter %d because rollback attempts ' ...
                    'reached the per-iter limit.\n'], iter);
                break;
            end

            oldGain = Gain;
            Gain = reduce_gain(Gain, config.Min_Gain);
            if abs(Gain - oldGain) < eps
                exitStatus = 'min_gain_rejected';
                fprintf(['  Stopping without accepting iter %d because Gain is already ' ...
                    'at the minimum %.3f.\n'], iter, config.Min_Gain);
                break;
            end
            fprintf('  Gain reduced from %.3f to %.3f.\n', oldGain, Gain);
            fprintf('  Recomputing Step4 for iter %d to overwrite correction file for iter %d.\n', ...
                iter - 1, iter);
            Step4_Iterative_Correction_Hetero(iter - 1, Gain);
            fprintf('  Rerunning iter %d with the overwritten correction file.\n', iter);
            continue;
        end

        rollbackAttemptsThisIter = 0;
        acceptedGain = Gain;
        currentScore = max([currentRMSE / config.Target_RMSE_K, ...
            currentMaxAbsError / config.Target_MaxAbsError_K]);
        fprintf('  Normalized convergence score = %.4f\n', currentScore);

        anyScoreImproved = currentScore < bestAccepted.score;
        meaningfulScoreImproved = isinf(bestAccepted.score) || ...
            currentScore < bestAccepted.score * (1 - config.Min_Score_Improvement);
        if anyScoreImproved
            bestAccepted.iter = iter;
            bestAccepted.rmse = currentRMSE;
            bestAccepted.maxAbsError = currentMaxAbsError;
            bestAccepted.score = currentScore;
            bestAccepted.gain = acceptedGain;
        end

        if meaningfulScoreImproved
            stalledAcceptedIterations = 0;
        else
            stalledAcceptedIterations = stalledAcceptedIterations + 1;
            fprintf('  No meaningful best-score improvement for %d accepted iter(s).\n', ...
                stalledAcceptedIterations);
        end

        hasConverged = ...
            currentRMSE <= config.Target_RMSE_K && currentMaxAbsError <= config.Target_MaxAbsError_K;

        if hasConverged
            exitStatus = 'converged';
            fprintf('\nConverged at iter %d.\n', iter);
            fprintf('  RMSE <= %.2f K and max local error is within approximately +/-%.2f K.\n', ...
                config.Target_RMSE_K, config.Target_MaxAbsError_K);
            lastAccepted.iter = iter;
            lastAccepted.rmse = currentRMSE;
            lastAccepted.maxAbsError = currentMaxAbsError;
            lastAccepted.score = currentScore;
            lastAccepted.gain = acceptedGain;
            break;
        end

        if currentRMSE > prevRMSE && abs(Gain - 0.8) < eps
            fprintf('\nGain shift: RMSE rebounded from %.4f K to %.4f K; next iter uses Gain = 0.5.\n', ...
                prevRMSE, currentRMSE);
            Gain = 0.5;
        end

        if maxPosError <= config.Target_MaxAbsError_K && abs(maxNegError) <= config.Target_MaxAbsError_K
            if Gain >= 0.5
                fprintf('\nGain shift: local errors are within +/-%.2f K; next iter uses fine-tuning Gain = 0.3.\n', ...
                    config.Target_MaxAbsError_K);
                Gain = 0.3;
            end
        end

        lastAccepted.iter = iter;
        lastAccepted.rmse = currentRMSE;
        lastAccepted.maxAbsError = currentMaxAbsError;
        lastAccepted.score = currentScore;
        lastAccepted.gain = acceptedGain;
        prevRMSE = currentRMSE;
        prevMaxAbsError = currentMaxAbsError;

        if stalledAcceptedIterations >= config.Stall_Patience
            exitStatus = 'stalled';
            fprintf('\nStopping because convergence score stalled for %d accepted iterations.\n', ...
                config.Stall_Patience);
            break;
        end

        iter = iter + 1;
    end

    if iter > config.Max_Iterations && strcmp(exitStatus, 'running')
        exitStatus = 'max_iterations';
        fprintf('\nStopped after reaching Max_Iterations = %d.\n', config.Max_Iterations);
    end

    fprintf('\nWorkflow finished.\n');
    fprintf('Exit status: %s\n', exitStatus);
    fprintf('Last accepted state: iter %d, RMSE %.4f K, max abs error %.4f K, Gain %.3f, score %.4f.\n', ...
        lastAccepted.iter, lastAccepted.rmse, lastAccepted.maxAbsError, ...
        lastAccepted.gain, lastAccepted.score);
    fprintf('Best accepted state: iter %d, RMSE %.4f K, max abs error %.4f K, Gain %.3f, score %.4f.\n', ...
        bestAccepted.iter, bestAccepted.rmse, bestAccepted.maxAbsError, ...
        bestAccepted.gain, bestAccepted.score);
end

function copy_resume_artifacts(sourceRunDir, runDir, startStep, startIter)
    startOrder = step_order(startStep);
    copied = {};

    if startOrder > step_order('Step1a')
        copied{end + 1} = copy_latest_required(sourceRunDir, runDir, ...
            'ICM_Hetero_YRow_*.mph', 'base Step1a COMSOL model');
    end

    if startOrder > step_order('Step1b')
        copied{end + 1} = copy_latest_required(sourceRunDir, runDir, ...
            'G_Matrix_Hetero_*.mat', 'Step1b G matrix');
    end

    if startOrder > step_order('Step2b')
        copied = [copied, copy_linewidth_input(sourceRunDir, runDir, startIter)];
    end

    if startOrder > step_order('Step3')
        copied{end + 1} = copy_latest_required(sourceRunDir, runDir, ...
            sprintf('Verify_Model_Hetero_*_diode_%d.mph', startIter), ...
            sprintf('Step3 verification model for iter %d', startIter));
        fprintf(['Step4 resume expects the copied verification model to contain ' ...
            'readable solution data.\n']);
    end

    fprintf('\nCopied resume artifacts into new run directory:\n');
    for k = 1:numel(copied)
        fprintf('  %s\n', copied{k});
    end
end

function copied = copy_linewidth_input(sourceRunDir, runDir, iter)
    if iter == 0
        copied = {copy_latest_required(sourceRunDir, runDir, ...
            'Initial_Optimization_Solution_Hetero_*.mat', ...
            'initial Step2b optimization solution')};
    else
        copied = {copy_latest_required(sourceRunDir, runDir, ...
            sprintf('Correction_Step_Latest_Hetero_*_diode_%d.mat', iter), ...
            sprintf('correction solution for iter %d', iter))};
    end
end

function destPath = copy_latest_required(sourceDir, destDir, pattern, description)
    matches = dir(fullfile(sourceDir, pattern));
    if isempty(matches)
        error('Cannot resume: missing %s matching %s in %s.', description, pattern, sourceDir);
    end
    [~, idx] = max([matches.datenum]);
    sourcePath = fullfile(matches(idx).folder, matches(idx).name);
    destPath = fullfile(destDir, matches(idx).name);
    copyfile(sourcePath, destPath);
end

function order = step_order(stepName)
    switch normalize_start_step(stepName)
        case 'Step1a'
            order = 1;
        case 'Step1b'
            order = 2;
        case 'Step2b'
            order = 3;
        case 'Step3'
            order = 4;
        case 'Step4'
            order = 5;
        otherwise
            error('Unknown step name: %s', stepName);
    end
end

function ensure_project_folders(projectRoot)
    folderNames = {'models', 'matlab', 'results', 'logs'};
    for k = 1:numel(folderNames)
        folderPath = fullfile(projectRoot, folderNames{k});
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
    end
end

function ensure_comsol_livelink_connected()
    mliPath = comsol_mli_path();
    serverPort = comsol_server_port();

    if exist(mliPath, 'dir') ~= 7
        error('COMSOL LiveLink mli path not found: %s', mliPath);
    end
    addpath(mliPath);

    mphstartStatus = exist('mphstart', 'file');
    if ~ismember(mphstartStatus, [2 6])
        error(['mphstart was not found after adding COMSOL LiveLink path: %s. ' ...
            'exist(''mphstart'', ''file'') returned %d.'], mliPath, mphstartStatus);
    end
    fprintf('Using mphstart from: %s\n', which('mphstart'));

    try
        mphstart(serverPort);
        fprintf('Connected to COMSOL Multiphysics Server on port %d.\n', serverPort);
    catch ME
        error(['Could not connect to COMSOL Multiphysics Server on port %d. ' ...
            'Start comsolmphserver.exe with -port %d before running this workflow. ' ...
            'Original error: %s'], serverPort, serverPort, ME.message);
    end
end

function [runDir, runName] = create_timestamped_run_dir(projectRoot, timestamp, prefix)
    baseRunName = [prefix '_' timestamp];
    runName = baseRunName;
    runDir = fullfile(projectRoot, 'results', runName);
    suffix = 1;

    while exist(runDir, 'dir')
        suffix = suffix + 1;
        runName = sprintf('%s_%02d', baseRunName, suffix);
        runDir = fullfile(projectRoot, 'results', runName);
    end

    mkdir(runDir);
end

function newGain = reduce_gain(currentGain, minGain)
    if currentGain > 0.5
        newGain = 0.5;
    elseif currentGain > 0.3
        newGain = 0.3;
    else
        newGain = max(currentGain * 0.5, minGain);
    end
end

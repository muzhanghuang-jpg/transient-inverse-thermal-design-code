function Run_Resume_Workflow(sourceRunDir, startStep, startIter)
%RUN_RESUME_WORKFLOW Resume the workflow in a new run directory.
%
% Examples:
%   Run_Resume_Workflow(fullfile(pwd, 'results', 'PRIOR_RUN'), 'Step2b')
%   Run_Resume_Workflow(fullfile(pwd, 'results', 'PRIOR_RUN'), 'Step3', 0)

    if nargin < 2 || isempty(startStep)
        startStep = 'Step2b';
    end
    if nargin < 3 || isempty(startIter)
        startIter = 0;
    end

    options = struct();
    options.Mode = 'resume';
    options.SourceRunDir = sourceRunDir;
    options.StartStep = startStep;
    options.StartIter = startIter;
    Run_Workflow_Core(options);

end

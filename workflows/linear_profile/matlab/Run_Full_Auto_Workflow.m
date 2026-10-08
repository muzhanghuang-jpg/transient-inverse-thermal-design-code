function Run_Full_Auto_Workflow()
%RUN_FULL_AUTO_WORKFLOW Timestamped COMSOL-MATLAB closed-loop pipeline.

    options = struct();
    options.Mode = 'full';
    Run_Workflow_Core(options);

end

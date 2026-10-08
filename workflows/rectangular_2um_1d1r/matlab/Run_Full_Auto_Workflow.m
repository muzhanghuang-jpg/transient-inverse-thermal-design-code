function Run_Full_Auto_Workflow()
    % --- 0. Material parameters ---
    k_SOI_val    = 130;  
    k_Hdoped_val = 1;    
    k_Ldoped_val = 1;    

    % --- 1. Response model and initial inverse solution ---
    % fprintf('\n[1/4] Build and solve the thermal response model.\n');
    % Step1a_Generate_ICM_Model_Hetero(k_SOI_val, k_Hdoped_val, k_Ldoped_val);
    % 
    % fprintf('\n[2/4] Extract the response matrix.\n');
    % Step1b_Extract_G_Matrices_6G_Decoupled(); 
    % 
    % fprintf('\n[3/4] Calculate the initial inverse solution (iteration 0).\n');
    % Step2b_optimization_Physical_WidthGen_Hetero();
    
    % --- 2. Nonlinear correction with gain reduction and rollback ---
    Max_Iterations = 10;
    Target_RMSE = 2.0; 
    
    current_Gain = 0.5; % Initial correction gain
    iter = 1;
    prev_max_error = inf;    % Maximum absolute error from the previous state
    prev_rmse = inf; 

    fprintf('\n[4/4] Start nonlinear verification and correction.\n');
    
    while iter <= Max_Iterations
        fprintf('\n==================================================\n');
        fprintf('Evaluate iteration %d with gain %.2f.\n', iter, current_Gain);
        fprintf('==================================================\n');
        
        % Build and solve the current layout.
        Step3_Build_Verify_Model_Hetero_diode(iter, k_SOI_val, k_Hdoped_val, k_Ldoped_val);
        
        % Evaluate the current field and prepare the next correction.
        [current_rmse, max_pos, max_neg] = Step4_Iterative_Correction_Hetero(iter, current_Gain);

        current_max_error= max(max_pos,abs(max_neg));
        
        fprintf('\nIteration %d results:\n', iter);
        fprintf('   - RMSE = %.2f K\n', current_rmse);
        fprintf('   - Maximum positive error = %+.2f K\n', max_pos);
        fprintf('   - Maximum negative error = %+.2f K\n', max_neg);
        
        % ---------------------------------------------------------
        % Reject an update when the maximum absolute error increases.
        % ---------------------------------------------------------
        if iter > 0 && current_max_error > prev_max_error 
            if current_Gain == 0.8
                fprintf('\nRollback: error increased with gain 0.8.\n');
                fprintf('Discard iteration %d and retry with gain 0.5.\n', iter);
                current_Gain = 0.5;
                Step4_Iterative_Correction_Hetero(iter - 1, current_Gain); % Recompute the update.
                continue; % Retry without advancing the iteration.
            elseif current_Gain == 0.5
                fprintf('\nRollback: error increased with gain 0.5.\n');
                fprintf('Discard iteration %d and retry with gain 0.3.\n', iter);
                current_Gain = 0.3;
                Step4_Iterative_Correction_Hetero(iter - 1, current_Gain); % Recompute with the reduced gain.
                continue; 
            end
        end
        
        % ---------------------------------------------------------
        % Stop when both sampled-error criteria are satisfied.
        % ---------------------------------------------------------
        if current_rmse <= Target_RMSE && max_pos < 5 && abs(max_neg) < 5
            fprintf('\nSampled-error criteria satisfied; stop correction.\n');
            break;
        end
        
        % ---------------------------------------------------------
        % Reduce the correction gain for the next update.
        % ---------------------------------------------------------
        % Check whether RMSE increased.
        if current_rmse > prev_rmse && current_Gain == 0.8
            fprintf('\nRMSE increased (%.2f -> %.2f); reduce gain to 0.5.\n', prev_rmse, current_rmse);
            current_Gain = 0.5;
        end
        
        % Reduce gain when both maximum signed errors are within 5 K.
        
        if max_pos <= 5 && abs(max_neg) <= 5
            if current_Gain >= 0.5
                fprintf('\nMaximum sampled errors are within 5 K; reduce gain to 0.3.\n');
                current_Gain = 0.3;
            
            end
        end
        
        % Advance to the next state.
        prev_max_error = current_max_error;
        prev_rmse = current_rmse;
        iter = iter + 1;
        
        if iter > Max_Iterations
            fprintf('\nMaximum iteration count reached; stop correction.\n');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the programmer: Abraham %
% Date: 2017-09-29                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose
% A generic implementation of Algorithm 1

%% Inputs
% xmax                   --- Maximizer of reachAvoidCostFunction
% direction              --- Direction to find the relative boundary
% thetamax               --- Upper bound of theta^\ast
% tolerance              --- Tolerance of theta^\ast (Exiting interval length)
% reachAvoidCostFunction --- Function handle to \hat{r}_{xmax + theta*direction}^U_vec( S, T)
% A_inequalities_input   --- U^N expressed as polytope (A,b)
% b_inequalities_input   --- U^N expressed as polytope (A,b)
% alpha_for_level_set    --- Threshold for the level set

%% Outputs
% thetaoptimal 


%% Notes
% Requires RAprob.m that serves as the oracle for computing 

function [optimal_theta, optimal_reachAvoid, optimal_inputs] = ...        
        getBoundaryPointViaBisectionNonFeas_CWH(xmax,...
                                                direction,...
                                                reachAvoidTube_A,...
                                                reachAvoidTube_b,...
                                                concatenated_A_matrix,...
                                                H_matrix_no_initial_state,...
                                                input_dimension,...
                                                last_time_step,...
                                                concatenated_state_sigma_no_input,...
                                                alpha_for_level_set,...
                                                tolerance_bisection,...
                                                lower_bound_on_optimal_theta,...
                                                upper_bound_on_optimal_theta,...
                                                initial_guess_for_U_at_phi_and_maxReachAvoid,...
                                                cost_function_of_phi_and_U,...
                                                PSoptions,...
                                                A_inequalities_input,...
                                                b_inequalities_input)
    optimal_reachAvoid_so_far = initial_guess_for_U_at_phi_and_maxReachAvoid(1);
    optimal_inputs_so_far = initial_guess_for_U_at_phi_and_maxReachAvoid(2:end);
    optimal_theta_so_far = 0;
    if (upper_bound_on_optimal_theta - lower_bound_on_optimal_theta) > tolerance_bisection
        fprintf('RA prob  |  Theta  | UstarEnergy | LB_theta | UB_theta | Exit reason (ExitFlag)\n');
        while (upper_bound_on_optimal_theta - lower_bound_on_optimal_theta) > tolerance_bisection
            initial_guess_for_input_vector = optimal_inputs_so_far;
            % Set phi as the middle point of the interval
            phi = (upper_bound_on_optimal_theta+lower_bound_on_optimal_theta)/2;
            % Slice constraint already enforced by direction and x_0 \in X by theta
            cost_function_at_phi = @(U) cost_function_of_phi_and_U([phi;U]);
            %%%%%%
            % Maximize \hat{r} subject to the control constraints
            % minimize -log(reachAvoidFunction(phi,U))
            % s.t.   input_vector \in U^N
            % at every phi for a U. 
            % If the max value is greater than alpha_for_level_set,
            % then feasible. Else, stop
            %%%%%%%%%%%%%%%%%%%%%%%%
            %PSoptions.Display = 'iter';
            [U_star,optimal_negative_log_reachAvoid_value,exitflag] = ...
                                                 patternsearch(cost_function_at_phi,...
                                                               initial_guess_for_input_vector,...
                                                               A_inequalities_input,...
                                                               b_inequalities_input,...
                                                               [],[],... % Enforces the slice
                                                               [],[],[],...
                                                               PSoptions);
            % Decode the exit message
            if exitflag > 0
                % Update the solution
                optimal_inputs_so_far = U_star;
                optimal_reachAvoid_so_far = exp(-optimal_negative_log_reachAvoid_value);
                %disp(output.message);                                             %Print message from patternsearch            
                if optimal_reachAvoid_so_far >= alpha_for_level_set
                    % Update theta to phi only if it is still within the "circle" (underapproximation set)
                    optimal_theta_so_far = phi;
                    exitmessage = '  Feasible';
                else
                    exitmessage = '  Infeasible';
                end            
            end
            % Print current solution
            fprintf(strcat('%1.4f  |  %1.4f  |  %1.4f  |  %1.4f  |  %1.4f  |  ',exitmessage,' \n'),...
                optimal_reachAvoid_so_far,...
                optimal_theta_so_far,...
                optimal_inputs_so_far'*optimal_inputs_so_far,...
                lower_bound_on_optimal_theta,...
                upper_bound_on_optimal_theta);
            % Update bounds
            if exitflag > 0 && optimal_reachAvoid_so_far >= alpha_for_level_set
                % Theta proved to be feasible
                lower_bound_on_optimal_theta = phi;
%                 upper_bound_on_optimal_theta = original_upper_bound;
            else
                % Shrink the upper bound
                upper_bound_on_optimal_theta = phi;
            end
        end        
    else
        fprintf('###Skipping bisection since theta interval too small! Returning the initial guess!\n')
    end
    optimal_inputs = optimal_inputs_so_far;
    optimal_reachAvoid = optimal_reachAvoid_so_far;
    optimal_theta = optimal_theta_so_far;
end

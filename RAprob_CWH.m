%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the programmer: Abraham %
% Date: 2018-03-23                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose
% Compute the reach-avoid probability using Genz's code

%% Inputs
% U                     : Input vector under evaluation
% concatenated_state_mean_without_initial_state : Abar * x_0 + G_matrix * disturbance_mean_concatenated_vector
% concatenated_state_sigma                      : Sigma of the concatenated state vector
% H_matrix_without_initial_state                : H_matrix * U + concatenated_state_mean_without_initial_state = concatenated_state_mean
% reachAvoidTube_A                              : Matrix A for the linear inequalities describing the polytopic reach-avoid tube
% reachAvoidTube_b                              : Vector b for the linear inequalities describing the polytopic reach-avoid tube
% state_dimension                               : n of the system
% myeps                                         : Tolerance for the Genz's algorithm's error estimate

%% Outputs
% prob                                          : \hat{r}_{\overline{x}_0}^U (Eqn 11/12 in the NAASS 2018 paper)

%% Notes
% Uses Genz's algorithm to compute the integral seen in (12)

function prob=RAprob_CWH(U,concatenated_state_mean_without_initial_state,concatenated_state_sigma,H_matrix_without_initial_state,reachAvoidTube_A,reachAvoidTube_b,state_dimension,myeps)    
    mu_vector=concatenated_state_mean_without_initial_state+H_matrix_without_initial_state*U;
    sigma_matrix=concatenated_state_sigma(state_dimension+1:end,state_dimension+1:end);
    %% QSCMVNV in a loop using the error estimate
    error_quadrature=10;
    points_base=10;
    points_power=1;
    while abs(error_quadrature)>myeps
        [prob,error_quadrature]=qscmvnv( points_base^points_power, sigma_matrix, repmat(-Inf,[size(reachAvoidTube_A,1),1]),reachAvoidTube_A, reachAvoidTube_b-reachAvoidTube_A*mu_vector);
        prob=round(prob/myeps)*myeps;
        % If zero, then saturate it for log
        if prob<myeps
            prob=myeps;
        end
        if points_power>5
            fprintf('Exceeded 5 iterations --- Required accuracy: %1.2e | Current error: %1.2e\n',myeps,error_quadrature);
%             break
        end        
        points_power=points_power+1;
    end
%    if points_power>5
%        fprintf('Took %d iterations\n',points_power-1)
%    end

%     %% Other methods
%     prob=qsimvnauto( sigma_matrix, reachAvoidTubeLB-mu_vector,  reachAvoidTubeUB-mu_vector,myeps,1e5);
%     prob=mvncdf(reachAvoidTubeLB,reachAvoidTubeUB,mu_vector,sigma_matrix);
%     prob=qscmvnv( 50000, sigma_matrix, reachAvoidTubeLB-mu_vector, eye(size(mu_vector,1),size(mu_vector,1)), reachAvoidTubeUB-mu_vector);
end

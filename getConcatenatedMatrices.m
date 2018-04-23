%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the programmer: Abraham %
% Date: 2018-03-23                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose
% For an arbitrary LTI dynamics, getConcatenatedMatrices produces the concatenated dynamics matrices.
% The matrices given by this function can be used to compute the concatenated state vector given a time horizon, control input vector, and the disturbance vector.

%% Inputs
% system_matrix     : A matrix
% input_matrix      : B matrix
% disturbance_matrix: F matrix
% last_time_step    : N 
%% Outputs
% concatenated_A_matrix : Matrix to be multiplied with initial state
% H_matrix              : Matrix to be multiplied with the input vector
% G_matrix              : Matrix to be multiplied with the disturbance vector

%% Notes
function [concatenated_A_matrix, H_matrix, G_matrix] = getConcatenatedMatrices(system_matrix, input_matrix, disturbance_matrix, last_time_step)
    % Dimensions of state, input, and disturbance
    state_dimension=size(system_matrix,2);
    input_dimension=size(input_matrix,2);
    disturbance_dimension=size(disturbance_matrix,2);

    %% concatenated_A matrix creation
    concatenated_A_matrix=eye(state_dimension);
    for rowNumber=1:last_time_step
        concatenated_A_matrix=[concatenated_A_matrix;system_matrix*concatenated_A_matrix(end-state_dimension+1:end,:)];
    end

    %% H creation
    H_matrix=zeros(state_dimension*(last_time_step+1),input_dimension*last_time_step);
    for rowNumber=1:last_time_step
        % Construct H_matrix block row wise --- n*(mT)
        H_matrix_temp=zeros(state_dimension,input_dimension*last_time_step);
        % What is the power of A in the first column? Ignoring the first block
        % row of zeros, it is the row number-1. (Ignoring is done later)
        maximum_exponent_for_system_matrix=rowNumber-1;
        % Construct the H_matrix block row of interest
        for system_matrix_exponent=maximum_exponent_for_system_matrix:-1:0
            column_left_indx=input_dimension*(maximum_exponent_for_system_matrix-system_matrix_exponent);
            H_matrix_temp(:,column_left_indx+1:column_left_indx+input_dimension)=system_matrix^system_matrix_exponent*input_matrix;
        end
        % The indices in the LHS ensures first row block is skipped
        H_matrix( rowNumber*state_dimension + 1 : rowNumber*state_dimension + state_dimension,:)=H_matrix_temp;    
    end

    %% G creation
    G_matrix=zeros(state_dimension*(last_time_step+1),disturbance_dimension*last_time_step);
    for rowNumber=1:last_time_step
        % Construct H_matrix block row wise --- n*(mT)
        G_matrix_temp=zeros(state_dimension,disturbance_dimension*last_time_step);
        % What is the power of A in the first column? Ignoring the first block
        % row of zeros, it is the row number-1. (Ignoring is done later)
        maximum_exponent_for_system_matrix=rowNumber-1;
        % Construct the H_matrix block row of interest
        for system_matrix_exponent=maximum_exponent_for_system_matrix:-1:0
            column_left_indx=disturbance_dimension*(maximum_exponent_for_system_matrix-system_matrix_exponent);
            G_matrix_temp(:,column_left_indx+1:column_left_indx+disturbance_dimension)=system_matrix^system_matrix_exponent*disturbance_matrix;
        end
        % The indices in the LHS ensures first row block is skipped
        G_matrix( rowNumber*state_dimension + 1 : rowNumber*state_dimension + state_dimension,:)=G_matrix_temp;    
    end
end

clear
clc
close all

params_CWH

%% Time convention
% Evolution starts from 0 and ends at time_horizon
% Hence, the total number of timesteps is time_horizon+1
% To make this clear, I will use last_time_step to refer the last time
% state is known 
last_time_step = time_horizon;

initial_xmax=initial_guess_for_xmax;
initial_guess_for_input_vector=zeros(last_time_step*input_dimension,1);

%% Initialization
polyarray_underapprox = [];
array_elapsed_time_PS = zeros(1,no_of_alphas);
optimal_theta_i_array = zeros(no_of_alphas,no_of_direction_vectors);
optimal_reachAvoid_i_array = zeros(no_of_alphas,no_of_direction_vectors);
optimal_inputs_i_cell = cell(no_of_alphas,1);

%% Compute concatenated matrices for state, input, and disturbance
%[x_0; x_1;..., x_T]=concatenated_A_matrix * x_0 + H_matrix * U + G_matrix * W;
[concatenated_A_matrix, H_matrix, G_matrix] = getConcatenatedMatrices(system_matrix, input_matrix, disturbance_matrix, last_time_step);
H_matrix_no_initial_state=H_matrix(state_dimension+1:end,:);

%% Create the reach-avoid open-loop cost function
% Reach-avoid tube constraint
% Stay in the safe set for t=1,...,T-1 and Reach the target set at T
% Safety constraint at t=0 is assumed to hold (TO BE ADDED IN CONSTRAINTS)
reachAvoidTube_A=blkdiag(kron(eye(last_time_step-1),A_safe_set),A_target_set);
reachAvoidTube_b=[kron(ones(last_time_step-1,1),b_safe_set);
                  b_target_set];
% U^N constraint
A_inequalities_input = kron(eye(last_time_step),A_input_set);
b_inequalities_input = kron(ones(last_time_step,1),b_input_set);
% Stochastics of G*W (initial state affects only the mean not the Sigma)
concatenated_disturbance_mean=kron(ones(last_time_step,1),mean_vector);
concatenated_disturbance_sigma=kron(eye(last_time_step),sigma_matrix);
concatenated_state_mean_no_input_and_no_initial_state_effect = G_matrix*concatenated_disturbance_mean;
concatenated_state_sigma_no_input=G_matrix*concatenated_disturbance_sigma*G_matrix';

timerval_polytopic_underapprox = tic;

%%% Guess for the patternsearch
%warning('off');
%cvx_begin quiet
%    variable initial_input_vector(input_dimension*last_time_step)
%    variable initial_xmax(state_dimension)
%    minimize norm(initial_input_vector)
%    subject to
%    reachAvoidTube_A*(concatenated_A_matrix(5:end,:)*initial_xmax + H_matrix_no_initial_state*initial_input_vector)<= reachAvoidTube_b 
%    -umax <= initial_input_vector <= umax 
%    A_safe_set*initial_xmax <= b_safe_set
%    [0 0 1 0; 0 0 0 1]*initial_xmax == slice_at_vx_vy' 
%cvx_end
%warning('on');
%initial_guess_for_xmax = initial_xmax;
%initial_guess_for_input_vector = initial_input_vector;

%% Compute the xmax using patternsearch
% Takes in [initial_state;input_vector] \in R^{state_dimension+input_dimension*last_time_step} for u_0,...u_{T-1}
% minimize -log(reachAvoidProb)
% s.t   input_vector \in U^N
cost_function_xmax = @(Y) -log(reachAvoidCostFunctionAssumingValidInitialState_CWH(Y(1:state_dimension),...
                                  Y(state_dimension+1:state_dimension+input_dimension*last_time_step),...
                                  concatenated_A_matrix,...
                                  concatenated_state_mean_no_input_and_no_initial_state_effect,...
                                  concatenated_state_sigma_no_input,...
                                  H_matrix_no_initial_state,...
                                  reachAvoidTube_A,...
                                  reachAvoidTube_b,...
                                  state_dimension,...
                                  myeps));

%% Patternsearch to compute xmax
% initial_guess_for_xmax defined in the beginning of main_CWH
initial_guess_for_initial_state_and_input_vector = [initial_guess_for_xmax;
                                                    initial_guess_for_input_vector];
A_inequalities_xmax = blkdiag(A_safe_set,A_inequalities_input);
b_inequalities_xmax = [b_safe_set;b_inequalities_input];
A_equalities_xmax = [0,0,1,0,zeros(1,input_dimension*last_time_step);
                     0,0,0,1,zeros(1,input_dimension*last_time_step)];
b_equalities_xmax = slice_at_vx_vy';
PSoptions = psoptimset('Display',display_string,'TolMesh',mesh_tolerance_for_patternsearch,'TolCon',constraint_tolerance_for_patternsearch);
[Y_star,optimal_negative_log_reachAvoid_value] = patternsearch(cost_function_xmax,...
                                                               initial_guess_for_initial_state_and_input_vector,...
                                                               A_inequalities_xmax,...
                                                               b_inequalities_xmax,...
                                                               A_equalities_xmax,b_equalities_xmax,... % Enforces the slice
                                                               [],[],[],...
                                                               PSoptions);
max_reachAvoid_value = exp(-optimal_negative_log_reachAvoid_value);
fprintf('Max reach-avoid probability: %1.2f\n',max_reachAvoid_value);
%% Setting up the variables
PSoptions = psoptimset('Display',display_string,'TolMesh',mesh_tolerance_for_patternsearch,'TolCon',constraint_tolerance_for_patternsearch,'MaxFunEvals',MaxFunEvals);
%% Compute the reach-avoid set using Algorithms 1 and 2 (HSCC 2018)
for alpha_index = 1: no_of_alphas
    optimal_theta_i = zeros(1, no_of_direction_vectors);
    optimal_reachAvoid_i = zeros(1, no_of_direction_vectors);
    optimal_inputs_i = zeros(input_dimension*last_time_step,no_of_direction_vectors);

    alpha_for_level_set=alpha_for_level_set_vector(alpha_index);
    fprintf('\nAnalyzing alpha: %1.2f (%d/%d)\n',alpha_for_level_set,alpha_index,no_of_alphas);

    if max_reachAvoid_value >= alpha_for_level_set
        xmax = Y_star(1:state_dimension);
        optimal_inputs_for_max_reachAvoid = Y_star(state_dimension+1:end);
        fprintf('Polytopic underapproximation exists\n');
            
        %% Bounds
        lower_bound_logRAP = log(alpha_for_level_set);
        upper_bound_logRAP = Inf;
        lower_bound_on_theta = 0;
        % Construct the first initial guess
        initial_guess_for_U_at_phi_and_maxReachAvoid = [max_reachAvoid_value;
                                                        optimal_inputs_for_max_reachAvoid];             
        %% Compute theta^ast_i for each d_i in set_of_direction_vectors
        for direction_index = 1: no_of_direction_vectors
            fprintf('Analyzing direction (shown transposed) :%d/%d\n',direction_index,no_of_direction_vectors);
            direction = set_of_directions_vectors(:,direction_index);
            disp(direction');
            A_times_direction = A_safe_set*direction;
            %% Compute theta_max for the given direction and update upper_bound_vector_theta_i(2)
            % theta_max is the solution to the following optimization problem
            % minimize -theta
            % s.t.      theta*(A_safe_set*direction) <= b_safe_set-A_safe_set*xmax
            cvx_begin quiet
                variable upper_bound_on_theta(1)
                minimize -upper_bound_on_theta
                subject to
                    upper_bound_on_theta*(A_safe_set*direction) <= b_safe_set-A_safe_set*xmax
            cvx_end
            fprintf('\bUpper bound of theta: %1.2f\n',upper_bound_on_theta);
    
            % equality constraint takes Y=[theta;input_vector]
            cost_function_of_phi_and_U = @(Y) -log(reachAvoidCostFunctionAssumingValidInitialState_CWH(xmax + Y(1)*direction,...
                                                  Y(2:1+input_dimension*last_time_step),...
                                                  concatenated_A_matrix,...
                                                  concatenated_state_mean_no_input_and_no_initial_state_effect,...
                                                  concatenated_state_sigma_no_input,...
                                                  H_matrix_no_initial_state,...
                                                  reachAvoidTube_A,...
                                                  reachAvoidTube_b,...
                                                  state_dimension,...
                                                  myeps));        
            [optimal_theta_i(direction_index), optimal_reachAvoid_i(direction_index), optimal_inputs_i(:,direction_index)] = ...
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
                                                    lower_bound_on_theta,...
                                                    upper_bound_on_theta,...
                                                    initial_guess_for_U_at_phi_and_maxReachAvoid,...
                                                    cost_function_of_phi_and_U,...
                                                    PSoptions,...
                                                    A_inequalities_input,...
                                                    b_inequalities_input);
        end
        vertices_of_polytopic_underapprox = xmax+optimal_theta_i.*set_of_directions_vectors;
        polytopic_underapprox = Polyhedron('V',[vertices_of_polytopic_underapprox]');
        minVRep(polytopic_underapprox);
        if size(polytopic_underapprox.V,1)~=(no_of_direction_vectors)
            fprintf('\n### Convex hull ate away few points!\n No. of points: %d\n',size(polytopic_underapprox.V,1));
        end
    else
        fprintf('No polytopic underapproximation exists\n');
        polytopic_underapprox = Polyhedron.emptySet(state_dimension);
    end
    %% Construct the polytope
    elapsed_time_PS = toc(timerval_polytopic_underapprox);
%     fprintf('Time taken: %1.4f\n', elapsed_time_PS);

    optimal_theta_i_array(alpha_index,:) = optimal_theta_i;
    optimal_reachAvoid_i_array(alpha_index,:) = optimal_reachAvoid_i;
    optimal_inputs_i_cell{alpha_index} = optimal_inputs_i;
    polyarray_underapprox = [polyarray_underapprox,polytopic_underapprox];
    array_elapsed_time_PS(alpha_index) = elapsed_time_PS;
end
% disp(array_elapsed_time_PS);
fprintf('This completes the polytopic computation.\n\n');

%% CDC 2017 (Gleason, Vinod, Oishi): Load the Lagrangian-based solution
disp('Lagrangian-based solution');
disp('=========================');
disp('You may obtain the codes to rerun this computation from ');
disp('https://hscl.unm.edu/wp-content/uploads/CDC17.zip.');
if code_flow_flag == 1
    load('CDC2017_solution_cwh_save.mat','DsetTemp');
    LagrangianSolution = DsetTemp.slice([3,4], slice_at_vx_vy);
    label_cells={'Safe Set','Target Set','Lagrangian','Chance-constrained','Algorithm 1'};    
elseif code_flow_flag == 2
    disp('### Lagrangian-based approach fails in this case due ');
    disp('to computational geometry issues.');
    LagrangianSolution = Polyhedron.emptySet(2);
    label_cells={'Safe Set','Chance-constrained','Algorithm 1'};    
end

%% CDC 2013 (Lesser, Oishi, Erwin): Load the chance-constrained-based solution
disp('');
disp('Chance-constrained based solution');
disp('=================================');
disp('You may rerun this computation via Lesser_cvx_chance_constrained.m');
if code_flow_flag == 1
    %% CDC 2017
    load('Lesser_CCC_Figure5.mat'); 
elseif code_flow_flag == 2
    %% CDC 2013 --- restricted to left half
    load('Lesser_CCC_Figure6.mat'); 
end

fprintf('\n\nComputation times are (Table 2)\n');
if code_flow_flag == 1
    figure(5);
elseif code_flow_flag == 2
    figure(6);
end
clf
hold on;
plot(polytope_safe_set.slice([3,4], slice_at_vx_vy), 'alpha', 1,'color', 'g')
if code_flow_flag == 1
    plot(polytope_target_set.slice([3,4], slice_at_vx_vy), 'color', 'k')
    plot(LagrangianSolution, 'color', 'y', 'LineStyle', '--','alpha',0.9)
    disp('Lagrangian');
    disp(14.5/60);      % From CDC 2017 paper
end
caxis([0 1])
colormap_matrix = colormap('copper');
color_step_size = 1/(size(colormap_matrix,1)-1);
colorpolytope = colormap_matrix(round((0.8-.3)/color_step_size)+1,:);
[~,donotplot_handle]=contourf(x01,x02,Prob-.3,[.8-.3 .8-.3]);
set(get(get(donotplot_handle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend 
plot(Polyhedron('lb',-ones(2,1),'ub',ones(2,1))*0.0001+[0;0.05], 'alpha',1,'color', colorpolytope,'LineStyle','-','LineWidth',0.001);
plot(projection(polyarray_underapprox,[1,2]),'alpha',0.7,'color','m');
set(gca,'Fontsize',20);
axis equal
xlabel('$x$','interpreter','latex','Fontsize',20);
ylabel('$y$','interpreter','latex','Fontsize',20);
leg=legend(label_cells);
if code_flow_flag == 2
    set(leg,'Location','BestOutside')
    x01 = x01(1):0.1:-0.7;
    axis([x01(1) x01(end) x02(1) x02(end)])
    set(gca,'YTick',x02(1:10:end))
end
box on;
grid on;
disp('CDC2013')
disp(timeSpent/60);
disp('Algorithm 1');
disp(array_elapsed_time_PS/60);

%vertex = polytopic_underapprox.V;
%hold on
%scatter(xmax(1),xmax(2),100,'ko','filled')
%scatter(vertices_of_polytopic_underapprox(1,:),vertices_of_polytopic_underapprox(2,:),100,'bo','filled')

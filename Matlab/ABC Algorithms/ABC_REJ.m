
%% ABC implementation 1 most basic - ABC REJECTION ALGORITHM:
% Population Monte Carlo simulation 
% based om the paper: 
% Aproximate Bayesian Computation usink Markov Chain Monte Carlo simulation: DREAM(ABC) 

% Once summary statistics have been defined we are left with finding those
% values of the parameter for which <= epsilon.


N  = 50;        % Number of different turin simulations.
Ns = 650;       % Number of time entries for each turin simulation. 
Bw = 4e9;       % Bandwidth (5Mhz).

%% --------- Generate "observed data" used as Y_obs -----------------------------

param_T       = 7.8e-9; 
param_G0      = 4.07e-9;  % db = -83.9);  
param_lambda  = 10e-12;   
sigma_N       = 1.673e-4;   %  same as sqrt(28e-9);

[P_Y_observed, t_observed] = sim_turin_matrix(N, Bw, Ns, param_T, param_G0, param_lambda, sigma_N);

% Do summary statistics on the observed data 
S_observed = sumStatMeanMoment(t_observed, P_Y_observed); 
% -------- END observed data  -------------------------------------------

%%
iterParamsUpdate = 50;
iterDistance = 1000;
accepted_distance = zeros(iterParamsUpdate,1);
accepted_params = zeros(4,iterDistance);
final_params = zeros(4,iterParamsUpdate*iterDistance);
epsilon = 0.005; % Acceptable error 
distance = epsilon + 1;
total_iter = 1;
%% --- Initial max/min conditions for parameters prior distribution ------
% a = min , b = max
 T_a = 7.8e-17; 
 T_b = 7.8e-3;  
 G0_a = 8.89e-15;    % Power gain (not in dB)
 G0_b = 4.07e-6;     % Power gain (not in dB)
 lambda_a = 10e-15;
 lambda_b = 10e-8;
 sigmaN_a = sqrt(28e-14); 
 sigmaN_b = sqrt(28e-4);
%% -----------------------------------------------------------------------
tic
%% ABC Rejection algorith with epsilon and limit update 
for params_updates = 1:iterParamsUpdate
    for i = 1:iterDistance
        while distance > epsilon 
            %% STEP 1: 
            % Sample parameter from a defined prior distribution:
            
        	% T (Reverberation time):
            param_T = T_a + (T_b-T_a).*rand(1,1); % generate one random number
        	% G0 (Reverberation gain)  
            param_G0 = (G0_a + (G0_b-G0_a).*rand(1,1)); % generate one random number within the given limits.
            % lambda ()  
            param_lambda = lambda_a + (lambda_b-lambda_a).*rand(1,1); % generate one random number within the given limits.
            % sigma_N (Variance noise floor)
            param_sigma_N = sigmaN_a + (sigmaN_b-sigmaN_a).*rand(1,1); % generate one random number within the given limits.
        
            %% STEP 2:
            % Simulate data using Turing model based on parameters from STEP 1 
            [P_Y, t] = sim_turin_matrix(N, Bw, Ns, param_T, param_G0, param_lambda, param_sigma_N);
        
            %% STEP 3:
            %  Do summary statistics on the simulated dataset: 
            S_simulated = sumStatMeanMoment(t, P_Y);
         
            %% STEP 4: calculate eucledian distance function 
            % Calculate the distance function (difference between simulated and observed summary statistics.     
            % Calculated based on the three summary statistics of 0th, 1st and 2nd moment.  
            % See formular on page ....
         
            SS1 = ((S_simulated(1) - S_observed(1))/S_simulated(4))^2;
            SS2 = ((S_simulated(2) - S_observed(2))/S_simulated(5))^2;
            SS3 = ((S_simulated(3) - S_observed(3))/S_simulated(6))^2;
         
            distance = ((SS1 + SS2 + SS3)^(0.5))/1e24; % a scaling problem here - needs investigation
        end  
        
        accepted_distance(i) = distance;
        distance = epsilon + 1; % "Reset" the distance to be higher than epsilon.
        % When a distance is withn acceptable limits (distance < epsilon)
        % continue and save the parameter
        accepted_params(1,i) =  param_T;
        accepted_params(2,i) =  param_G0;
        accepted_params(3,i) =  param_lambda;
        accepted_params(4,i) =  param_sigma_N;
        
        % Save the accepted final parameter from each iteration
        final_params(1,total_iter) =  param_T;
        final_params(2,total_iter) =  param_G0;
        final_params(3,total_iter) =  param_lambda;
        final_params(4,total_iter) =  param_sigma_N;
        
        total_iter = total_iter + 1;
    end 
   
    % Update the parameter limits from the max/min values from the above
    % accepted parameters
    T_a      = min(accepted_params(1,:));
    T_b      = max(accepted_params(1,:));
    G0_a     = min(accepted_params(2,:));
    G0_b     = max(accepted_params(2,:));
    lambda_a = min(accepted_params(3,:));
    lambda_b = max(accepted_params(3,:));
    sigmaN_a = min(accepted_params(4,:));
    sigmaN_b = max(accepted_params(4,:));
    
    % Update epsilon: Need a better solution based on the research on the
    % epsilon update
    if epsilon > 0 
        epsilon = epsilon - 0.00001;
    end
    
    disp(params_updates)
end    
toc

x = 1:61;
x = x * 20;
t = tiledlayout(2,2, 'TileSpacing', 'none', 'Padding', 'compact');
title(t,'Parameter estimation using ABC rejection algorithm for Turin model');

nexttile
scatter(x,final_params(1,:))
yline(7.8e-9); % The actual parameter
title('T - reverbation time ')
% xtitle('Sample number');

nexttile 
scatter(x,final_params(2,:))
yline(4.07e-9);
title('G_0 - reverbation gain')
% xtitle('Sample number');


nexttile
scatter(x,final_params(3,:))
yline(10e-12);
title('\lambda - arrival rate')
% xtitle('Sample number');

nexttile
scatter(x,final_params(4,:))
yline(1.673e-4);
title('\sigma_n - sigma noise')
% xtitle('Sample number');



 
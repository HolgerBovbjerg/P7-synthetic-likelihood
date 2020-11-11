%% ABC implementation 2 - ABC REJECTION ALGORITHM:
clear
%% --- Global turin model simulation parameters ---------------------------------------------------
N  = 50;    % Number of different turin simulations.
Ns = 801;   % Number of time entries for each turin simulation. 
Bw = 4e9;   % Bandwidth (4Ghz).

%% --- Generate "observed data" -----------------------------------------------------

load("Theta_true_values.mat") 
% If Theta_true_values.mat not generated - uncoment and run the following block:

% T       = 7.8e-9;    
% G0      = db2pow(-83.9);
% lambda  = 10e9;
% sigma_N = 1.673e-4;
% M = 2000; % Number of summary statisctics realisations

% Theta_true_values = [T G0 lambda sigma_N];

S_obs = zeros(2000,4);

parfor i = 1:2000
    [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, Theta_true_values);
    S_obs(i,:) = create_statistics(Pv, t);
end
%%
mu_S_obs = mean(S_obs);     % Mean of the summary statistics 
Sigma_S_obs = cov(S_obs);     % Covariance of summary statistics

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
% a = min , b = max
 T_min = 1e-9; 
 T_max = 15e-9;  
 G0_min = db2pow(pow2db(theta_true(2)) - 10);    % Power gain (not in dB)
 G0_max = db2pow(pow2db(theta_true(2)) + 10);     % Power gain (not in dB)
 lambda_min = 1e8;
 lambda_max = 20e9;
 sigmaN_min = sqrt(0.28e-10); 
 sigmaN_max = sqrt(0.28e-8);
 
%% --- ABC rejection algorithm ---------------------------------------------------------------------
% Set total iterations
iterations = 10;

% Number of summary statistics sets to generate  
sumstat_iter = 2000;

% Extract this amount of parameter entries from each generated summary
% statistic
nbr_extract = 100;

params_T = zeros(iterations,nbr_extract);
params_G0 = zeros(iterations,nbr_extract);
params_lambda = zeros(iterations,nbr_extract);
params_sigma_N = zeros(iterations,nbr_extract);

disp('ABC algorithm computing, please wait... ')
tic

% Iteration 1
out = zeros(5,sumstat_iter);
d = zeros(sumstat_iter,1);
parfor i = 1:sumstat_iter
    %% STEP 1: Sample parameter from predefined prior distribution (uniform):      
    % T (Reverberation time):
    param_T = T_min + (T_max-T_min)*rand; % generate one random number
    % G0 (Reverberation gain)  
    param_G0 = (G0_min + (G0_max-G0_min)*rand); % generate one random number within the given limits.
    % lambda ()  
    param_lambda = lambda_min + (lambda_max-lambda_min)*rand; % generate one random number within the given limits.
    % sigma_N (Variance noise floor)
    param_sigma_N = sigmaN_min + (sigmaN_max-sigmaN_min)*rand; % generate one random number within the given limits.

    theta_curr = [param_T param_G0 param_lambda param_sigma_N];
    
    %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
    [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_curr);
    S_simulated = create_statistics(Pv, t);
    %% STEP 3: calculate the difference between observed and simulated summary statistics 
    % Mahalanobis distance see formular in document.
    d(i) = (S_simulated - mu_S_obs)/Sigma_S_obs * (S_simulated - mu_S_obs)';

    % Row 1 of the out vector contains the distance 
    % the rest of the rows contains the corresponding parameters 
    % used for generating that specific distance.    
    out(:,i) =  [d(i);...
                param_T;...
                param_G0;...
                param_lambda;...
                param_sigma_N];
    disp(i);
end 
% Sort the "out" matrix so that the lowest euclidean distance is at the
% (1,1) matrix position and highest distance is at (max,1) 
out = sortrows(out',1)';
%%
% Following vectors holds ALL extracted parameter values 
% that was within euclidean distance
params_T(1,:)       = out(2,1:nbr_extract);
params_G0(1,:)      = out(3,1:nbr_extract);
params_lambda(1,:)  = out(4,1:nbr_extract);
params_sigma_N(1,:) = out(5,1:nbr_extract);

% Update the prior for the next iteration
% [f_T,xi_T] = ksdensity(params_T(1,:));
% [f_G0,xi_G0] = ksdensity(params_G0(1,:));
% [f_lambda,xi_lambda] = ksdensity(params_lambda(1,:));
% [f_sigma_N,xi_sigma_N] = ksdensity(params_sigma_N(1,:));

for a = 2:iterations 
    out = zeros(5,sumstat_iter);
    d = zeros(sumstat_iter,1);
    parfor i = 1:sumstat_iter
        %% STEP 1: Sample parameter from predefined prior distribution (uniform):      
%         % T (Reverberation time):
%         param_T = randpdf(f_T,xi_T,[1, 1]); % generate one random number
%         % G0 (Reverberation gain)  
%         param_G0 = randpdf(f_G0,xi_G0,[1, 1]); % generate one random number within the given limits.
%         % lambda ()  
%         param_lambda = randpdf(f_lambda,xi_lambda,[1, 1]); % generate one random number within the given limits.
%         % sigma_N (Variance noise floor)
%         param_sigma_N = randpdf(f_sigma_N,xi_sigma_N,[1, 1]); % generate one random number within the given limits.
         
        param_T = T_min + (T_max-T_min)*rand; % generate one random number
        % G0 (Reverberation gain)  
        param_G0 = (G0_min + (G0_max-G0_min)*rand); % generate one random number within the given limits.
        % lambda ()  
        param_lambda = lambda_min + (lambda_max-lambda_min)*rand; % generate one random number within the given limits.
        % sigma_N (Variance noise floor)
        param_sigma_N = sigmaN_min + (sigmaN_max-sigmaN_min)*rand; % generate one random number within the given limits.
        
%         if param_T < T_min
%             param_T = T_min;
%         end
%         if param_T > T_max
%             param_T = T_max;
%         end
%         if param_G0 < G0_min
%             param_G0 = G0_min;
%         end
%         if param_T > G0_max
%             param_G0 = G0_max;
%         end
%         if param_lambda < lambda_min
%             param_lambda = lambda_min;
%         end
%         if param_lambda > lambda_max
%             param_lambda = lambda_max;
%         end
%         if param_sigma_N < sigmaN_min
%             param_sigma_N = sigmaN_min;
%         end
%         if param_sigma_N > sigmaN_max
%             param_sigma_N = sigmaN_max;
%         end
        
        theta_curr = [param_T param_G0 param_lambda param_sigma_N];
        %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
        
        [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_curr);
    	S_simulated = create_statistics(Pv, t);
        %% STEP 3: calculate the difference between observed and simulated summary statistics 
        % Mahalanobis distance
        d(i) = (S_simulated - mu_S_obs)/Sigma_S_obs * (S_simulated - mu_S_obs)';
        
        % Row 1 of the out vector contains the distance 
        % the rest of the rows contains the corresponding parameters 
        % used for generating that specific distance.    
        out(:,i) =  [d(i);...
                    param_T;...
                    param_G0;...
                    param_lambda;...
                    param_sigma_N];
        disp(i);
    end 
    % Sort the "out" matrix so that the lowest distance is at the
    % (1,1) matrix position and highest distance is at (1,end) 
    out = sortrows(out',1)';

    % Following vectors holds ALL extracted parameter values 
    % that was within euclidean distance
    params_T(a,:)       = out(2,1:nbr_extract);
    params_G0(a,:)      = out(3,1:nbr_extract);
    params_lambda(a,:)  = out(4,1:nbr_extract);
    params_sigma_N(a,:) = out(5,1:nbr_extract);
    
%    % Update the prior for the next iteration
%    [f_T,xi_T] = ksdensity(params_T(a,:));
%    [f_G0,xi_G0] = ksdensity(params_G0(a,:));
%    [f_lambda,xi_lambda] = ksdensity(params_lambda(a,:));
%    [f_sigma_N,xi_sigma_N] = ksdensity(params_sigma_N(a,:));
   
   T_min      = min(params_T(a,:));
   T_max      = max(params_T(a,:));
   G0_min     = min(params_G0(a,:));
   G0_max     = max(params_G0(a,:));
   lambda_min = min(params_lambda(a,:));
   lambda_max = max(params_lambda(a,:));
   sigmaN_min = min(params_sigma_N(a,:));
   sigmaN_max = max(params_sigma_N(a,:));
   disp(a);
end 
toc
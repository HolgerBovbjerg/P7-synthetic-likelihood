%% ABC implementation 2 - ABC REJECTION ALGORITHM:
clear
%% --- Global turin model simulation parameters ---------------------------------------------------
N  = 300;    % Number of different turin simulations.
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

load('S_obs.mat')
% load('S_obs_meas_data.mat')
%
% S_obs = zeros(2000,9);
%
% parfor i = 1:2000
%     [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_true);
%      S_obs(i,:) = create_statistics(Pv, t);
% end
%%
mu_S_obs = mean(S_obs);     % Mean of the summary statistics
Sigma_S_obs = cov(S_obs);     % Covariance of summary statistics

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
load('Prior_data_large_prior_min_max_values.mat')

%% --- ABC rejection algorithm ---------------------------------------------------------------------
% Set total iterations
iterations = 3;

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
    param_T = prior(1,1) + (prior(1,2) - prior(1,1)).*rand; % generate one random number
    % G0 (Reverberation gain)
    param_G0 = prior(2,1) + (prior(2,2) - prior(2,1)).*rand; % generate one random number within the given limits.
    % lambda ()
    param_lambda = prior(3,1) + (prior(3,2) - prior(3,1)).*rand; % generate one random number within the given limits.
    % sigma_N (Variance noise floor)
    param_sigma_N = prior(4,1) + (prior(4,2) - prior(4,1)).*rand; % generate one random number within the given limits.
    
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
%
% Following vectors holds ALL extracted parameter values
% that was within euclidean distance
params_T(1,:)       = out(2,1:nbr_extract);
params_G0(1,:)      = out(3,1:nbr_extract);
params_lambda(1,:)  = out(4,1:nbr_extract);
params_sigma_N(1,:) = out(5,1:nbr_extract);

T_min = min(params_T(1,:));
if T_min < prior(1,1)
    T_min = prior(1,1);
end
G0_min = min(params_G0(1,:));
if G0_min < prior(2,1)
    G0_min = prior(2,1);
end
lambda_min = min(params_lambda(1,:));
if lambda_min < prior(3,1)
    lambda_min = prior(3,1);
end
sigma_N_min = min(params_sigma_N(1,:));
if sigma_N_min < prior(4,1)
    sigma_N_min = prior(4,1);
end
T_max = max(params_T(1,:));
if T_max < prior(1,2)
    T_max = prior(1,2);
end
G0_max = max(params_G0(1,:));
if G0_max < prior(2,2)
    G0_max = prior(2,2);
end
lambda_max = max(params_lambda(1,:));
if lambda_max < prior(3,2)
    lambad_max = prior(3,2);
end
sigma_N_max = max(params_sigma_N(1,:));
if sigma_N_max < prior(4,2)
    sigma_N_max = prior(4,2);
end

%% Iterate

for a = 2:iterations
    out = zeros(5,sumstat_iter);
    d = zeros(sumstat_iter,1);
    parfor i = 1:sumstat_iter
        %% STEP 1: Sample parameter from predefined prior distribution (uniform):
        % T (Reverberation time):
        param_T = T_min + (T_max - T_min).*rand; % generate one random number
        % G0 (Reverberation gain)
        param_G0 = G0_min + (G0_max - G0_min).*rand; % generate one random number within the given limits.
        % lambda ()
        param_lambda = lambda_min + (lambda_max - lambda_min).*rand; % generate one random number within the given limits.
        % sigma_N (Variance noise floor)
        param_sigma_N = sigma_N_min + (sigma_N_max - sigma_N_min).*rand; % generate one random number within the given limits.
        
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
    
    T_min = min(params_T(a,:));
    if T_min < prior(1,1)
        T_min = prior(1,1);
    end
    G0_min = min(params_G0(a,:));
    if G0_min < prior(2,1)
        G0_min = prior(2,1);
    end
    lambda_min = min(params_lambda(a,:));
    if lambda_min < prior(3,1)
        lambda_min = prior(3,1);
    end
    sigma_N_min = min(params_sigma_N(a,:));
    if sigma_N_min < prior(4,1)
        sigma_N_min = prior(4,1);
    end
    T_max = max(params_T(a,:));
    if T_max > prior(1,2)
        T_max = prior(1,2);
    end
    G0_max = max(params_G0(a,:));
    if G0_max > prior(2,2)
        G0_max = prior(2,2);
    end
    lambda_max = max(params_lambda(a,:));
    if lambda_max > prior(3,2)
        lambad_max = prior(3,2);
    end
    sigma_N_max = max(params_sigma_N(a,:));
    if sigma_N_max > prior(4,2)
        sigma_N_max = prior(4,2);
    end
    disp(a);
end
disp('ABC algorithm finished... ')
toc
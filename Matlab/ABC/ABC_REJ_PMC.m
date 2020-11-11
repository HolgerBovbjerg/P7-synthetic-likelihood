%% ABC implementation 2 - ABC REJECTION ALGORITHM:
clear
%% --- Global turin model simulation parameters ---------------------------------------------------
N  = 50;    % Number of different turin simulations.
Ns = 801;   % Number of time entries for each turin simulation. 
Bw = 4e9;   % Bandwidth (4Ghz).

%% --- Generate "observed data" -----------------------------------------------------
param_T_obs       = 7.8e-9;  % Reverberation time
param_G0_obs      = db2pow(-83.9);    % linear gain 
param_lambda_obs  = 10e9;       % arrival rate (1/s)    
param_sigma_N_obs       = sqrt(0.28e-9);   % Noise std
M = 200; %number of summary statisctics realisations

% S_obs = create_statistics(M, N, param_T_obs , param_G0_obs,...
% param_lambda_obs, param_sigma_N_obs, Bw, Ns); % create new summary statistic of
% observed data
load S_obs;
%%
mu_S_obs = mean(S_obs);
Sigma_S_obs = cov(S_obs);

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
% a = min , b = max

G0_a = db2pow(param_G0_obs-10);    % Power gain (not in dB)
G0_b = db2pow(param_G0_obs+10);     % Power gain (not in dB)
tmax = (Ns-1)/Bw;
T_a = 0; 
T_b = (1/2)*tmax/100;  
lambda_a = 1/tmax;
lambda_b = 20e9;
sigmaN_a = sqrt(0.28e-8); 
sigmaN_b = sqrt(0.28e-10);
 
%% --- ABC rejection algorithm ---------------------------------------------------------------------
% Set total iterations
iterations = 3;

% Number of summary statistics sets to generate  
sumstat_iter = 20;

% Extract this amount of parameter entries from each generated summary
% statistic
nbr_extract = 10;

% Preallocation of vectors: 

params_T = zeros(iterations,nbr_extract);
params_G0 = zeros(iterations,nbr_extract);
params_lambda = zeros(iterations,nbr_extract);
params_sigma_N = zeros(iterations,nbr_extract);

param_T= zeros(sumstat_iter,1);
param_G0 = zeros(sumstat_iter,1);
param_lambda = zeros(sumstat_iter,1);
param_sigmaN= zeros(sumstat_iter,1);

disp('ABC algorithm computing, please wait... ')
tic

%% STEP 1: Sample parameters from predefined prior distribution (uniform):      
a = 3;
out = zeros(5,sumstat_iter);
d = zeros(sumstat_iter,1);
for i = 1:sumstat_iter
    %% STEP 1: Sample parameter from predefined prior distribution (uniform):      
    % T (Reverberation time):
    param_T(i) = T_a + (T_b-T_a)*rand; % generate one random number
    % G0 (Reverberation gain)  
    param_G0(i) = (G0_a + (G0_b-G0_a)*rand); % generate one random number within the given limits.
    % lambda ()  
    param_lambda(i) = lambda_a + (lambda_b-lambda_a)*rand; % generate one random number within the given limits.
    % sigma_N (Variance noise floor)
    param_sigma_N(i) = sigmaN_a + (sigmaN_b-sigmaN_a)*rand; % generate one random number within the given limits.

    %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
    S_simulated = create_statistics(1, N, param_T , param_G0, param_lambda, param_sigma_N, Bw, Ns);
    %% STEP 3: calculate the difference between observed and simulated summary statistics 
    % Mahalanobis distance see formular in document.
    d(i) = (S_simulated - mu_S_obs)/Sigma_S_obs * (S_simulated - mu_S_obs)';

    % Row 1 of the out vector contains the distance 
    % the rest of the rows contains the corresponding parameters 
    % used for generating that specific distance.    
    out(:,i) =  [d(i);...
                param_T(i);...
                param_G0(i);...
                param_lambda(i);...
                param_sigma_N(i)];
    disp(i);
end 
% Sort the "out" matrix so that the lowest euclidean distance is at the
% (1,1) matrix position and highest distance is at (max,1) 
out = sortrows(out',1)';

params = out(2:end,1:nbr_extract);
params_T(a,:)       = out(2,1:nbr_extract);
params_G0(a,:)      = out(3,1:nbr_extract);
params_lambda(a,:)  = out(4,1:nbr_extract);
params_sigma_N(a,:) = out(5,1:nbr_extract);


% Weights and variances
weights = ones(iterations, nbr_extract)./nbr_extract;
Var_T      = var(params_T(a,:));
Var_G0     = var(params_G0(a,:));
Var_lambda = var(params_lambda(a,:));
Var_sigmaN = var(params_sigma_N(a,:));

%% Update prior and iterate
param_T_star = zeros(sumstat_iter,1);
param_G0_star = zeros(sumstat_iter,1);
param_lambda_star = zeros(sumstat_iter,1);
param_sigmaN_star = zeros(sumstat_iter,1);
for a = 2:iterations 
    out = zeros(5,sumstat_iter);
    d = zeros(sumstat_iter,1);
    % T (Reverberation time):
    for j = 1:sumstat_iter
        param_T_star(j) = params_T(a-1,randi(nbr_extract)); % generate one random number
        % G0 (Reverberation gain)  
        param_G0_star(j) = params_G0(a-1,randi(nbr_extract)); % generate one random number within the given limits.
        % lambda ()  
        param_lambda_star(j) = params_lambda(a-1,randi(nbr_extract)); % generate one random number within the given limits.
        % sigma_N (Variance noise floor)
        param_sigmaN_star(j) = params_sigma_N(a-1,randi(nbr_extract)); % generate one random number within the given limits.
    end
    for i = 1:sumstat_iter
        %% STEP 1: Sample parameter from predefined prior distribution (uniform):     
        
        param_T = normrnd(param_T_star(i),Var_T);
        param_G0 = normrnd(param_G0_star(i),Var_G0);
        param_lambda = normrnd(param_lambda_star(i),Var_lambda);
        param_sigmaN = normrnd(param_sigmaN_star(i),Var_sigmaN);
        
        %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
        S_simulated = create_statistics(1, N, param_T , param_G0, param_lambda, param_sigma_N, Bw, Ns);
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
    end 
    
    % Sort the "out" matrix so that the lowest euclidean distance is at the
    % (1,1) matrix position and highest distance is at (max,1) 
    out = sortrows(out',1)';

    params = out(2:end,1:nbr_extract);
    params_T(a,:)       = out(2,1:nbr_extract);
    params_G0(a,:)      = out(3,1:nbr_extract);
    params_lambda(a,:)  = out(4,1:nbr_extract);
    params_sigma_N(a,:) = out(5,1:nbr_extract);
    
    Var_T      = var(params_T(a,:));
    Var_G0     = var(params_G0(a,:));
    Var_lambda = var(params_lambda(a,:));
    Var_sigmaN = var(params_sigma_N(a,:));
end 
toc
   



 
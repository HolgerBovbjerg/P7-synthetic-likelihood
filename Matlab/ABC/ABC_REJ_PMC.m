%% ABC-PMC ALGORITHM:
% As described in:
% Bharti, A., & Pedersen, T. (2020). Calibration of Stochastic Channel Models using Approximate Bayesian Computation.
% https://doi.org/10.1109/GCWkshps45667.2019.9024563
% (without regression)
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
% theta_true = [T G0 lambda sigma_N];

% M = 2000; % Number of summary statisctics realisations

% Theta_true_values = [T G0 lambda sigma_N];

% 
% [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_true);
% S_obs(i,:) = create_statistics(Pv, t);


load('S_obs.mat')
%%
mu_S_obs = mean(S_obs);     % Mean of the summary statistics
Sigma_S_obs = cov(S_obs);     % Covariance of summary statistics

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
load('Prior_data_large_prior_min_max_values.mat')

%% --- BSL PMC algorithm ---------------------------------------------------------------------
for i = 1:L
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr);
    s_sim(i,:) = create_statistics(Pv, t);
end

loglikelihood = synth_loglikelihood(s_obs,s_sim);

% Probability factor for generating population pool
prob_factor = 1e4;

params_T = zeros(iterations,nbr_extract);
params_G0 = zeros(iterations,nbr_extract);
params_lambda = zeros(iterations,nbr_extract);
params_sigma_N = zeros(iterations,nbr_extract);

disp('BSL algorithm computing, please wait... ')
tic

%% Iteration 1 - rejection step on uniform prior
out = zeros(5,sumstat_iter);
d = zeros(sumstat_iter,1);

% STEP 1: Sample parameter from predefined prior distribution (uniform):
% T (Reverberation time):
param_T = prior(1,1) + (prior(1,2) - prior(1,1)).*rand(sumstat_iter,1); % generate one random number
% G0 (Reverberation gain)
param_G0 = prior(2,1) + (prior(2,2) - prior(2,1)).*rand(sumstat_iter,1); % generate one random number within the given limits.
% lambda ()
param_lambda = prior(3,1) + (prior(3,2) - prior(3,1)).*rand(sumstat_iter,1); % generate one random number within the given limits.
% sigma_N (Variance noise floor)
param_sigma_N = prior(4,1) + (prior(4,2) - prior(4,1)).*rand(sumstat_iter,1); % generate one random number within the given limits.

parfor i = 1:sumstat_iter
    theta_curr = [param_T(i) param_G0(i) param_lambda(i) param_sigma_N(i)];
    
    % STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
    [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_curr);
    S_simulated = create_statistics(Pv, t);
    % STEP 3: calculate the difference between observed and simulated summary statistics
    % Mahalanobis distance see formular in document.
    d(i) = (S_simulated - mu_S_obs)/Sigma_S_obs * (S_simulated - mu_S_obs)';
    
    % Row 1 of the out vector contains the distance
    % the rest of the rows contains the corresponding parameters
    % used for generating that specific distance.
    out(:,i) =  [d(i);...
                theta_curr'];
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
accepted_params = [params_T; params_G0; params_lambda; params_sigma_N];

% Calculate first weights and covariance
weights_T = ones(1,nbr_extract)./nbr_extract;
weights_G0 = ones(1,nbr_extract)./nbr_extract;
weights_lambda = ones(1,nbr_extract)./nbr_extract;
weights_sigma_N = ones(1,nbr_extract)./nbr_extract;
weights = [weights_T; weights_G0; weights_lambda; weights_sigma_N];
var_T = var(params_T(1,:));
var_G0 = var(params_G0(1,:));
var_lambda = var(params_lambda(1,:));
var_sigma_N = var(params_sigma_N(1,:));
covariance = 2*diag([var_T var_G0 var_lambda var_sigma_N]);

% Choose theta from accepted parameters of last iteration with propbability based on wieghts
index_T = randsample((1:nbr_extract),sumstat_iter,true,weights(1,:));
index_G0 = randsample((1:nbr_extract),sumstat_iter,true,weights(2,:));
index_lambda = randsample((1:nbr_extract),sumstat_iter,true,weights(3,:));
index_sigma_N = randsample((1:nbr_extract),sumstat_iter,true,weights(4,:));

theta_prop = [params_T(1,index_T); params_G0(1,index_G0); params_lambda(1,index_lambda); params_sigma_N(1,index_sigma_N)];

%% sequential BSL Iterations (PMC)
for a = 1:20%iterations
    out = zeros(5,sumstat_iter);
    d = zeros(sumstat_iter,1);   
    parfor i = 1:sumstat_iter
        % Perturb theta 
        theta_curr = mvnrnd(theta_prop(:,i),covariance);
        while(check_params(theta_curr,prior)==2)
            theta_curr = mvnrnd(theta_prop(:,i),covariance);
        end
        theta_curr(3) = round(theta_curr(3));
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
                    theta_curr'];
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

    accepted_params = [params_T(a,:); params_G0(a,:); params_lambda(a,:); params_sigma_N(a,:)];
    
    % Compute weights for next iteration
    old_weights = weights; % Save current weights
    for l = 1:size(weights,1)
        for k = 1:size(weights,2)
            weights(l,k) = 1 /(sum(pdf('normal',accepted_params(l,:), accepted_params(l,k), sqrt(covariance(l,l)).* old_weights(l,:))));    
        end
    end
    weights = weights./sum(weights,2);
    % Find oovariance for next iteration
    covariance = 2*diag(diag(cov(out(2:5,1:nbr_extract)')));
    % Generate vector with each entry proportional weight for sampling 
    probs = round(weights*prob_factor); 
    % Generate population pool for sampling of new parameters
    big_T = [];
    big_G0 = [];
    big_lambda = [];
    big_sigma_N = [];
    for j = 1:nbr_extract
        big_T = [big_T repelem(accepted_params(1,j), probs(1,j))];
        big_G0 = [big_G0 repelem(accepted_params(2,j), probs(2,j))];
        big_lambda = [big_lambda repelem(accepted_params(3,j), probs(3,j))];
        big_sigma_N = [big_sigma_N repelem(accepted_params(4,j), probs(4,j))];
    end
    
    % Sample new para,eter values from population pool
    theta_prop = [datasample(big_T, sumstat_iter);...
                 datasample(big_G0, sumstat_iter);...
                 datasample(big_lambda, sumstat_iter);...
                 datasample(big_sigma_N, sumstat_iter)];

    disp(a); % display iteration number
end
disp('BSL algorithm finished... ')
toc
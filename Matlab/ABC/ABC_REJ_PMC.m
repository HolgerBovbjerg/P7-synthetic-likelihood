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
% theta_true = [T G0 lambda sigma_N];

% M = 2000; % Number of summary statisctics realisations

% Theta_true_values = [T G0 lambda sigma_N];

% S_obs = zeros(2000,9);
%
% parfor i = 1:2000
%     [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_true);
%      S_obs(i,:) = create_statistics(Pv, t);
% end

load('S_obs.mat')
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
param_T = zeros(sumstat_iter,1);
param_G0 = zeros(sumstat_iter,1);
param_lambda = zeros(sumstat_iter,1);
param_sigma_N = zeros(sumstat_iter,1);

parfor i = 1:sumstat_iter
    % STEP 1: Sample parameter from predefined prior distribution (uniform):
    % T (Reverberation time):
    param_T(i) = prior(1,1) + (prior(1,2) - prior(1,1)).*rand; % generate one random number
    % G0 (Reverberation gain)
    param_G0(i) = prior(2,1) + (prior(2,2) - prior(2,1)).*rand; % generate one random number within the given limits.
    % lambda ()
    param_lambda(i) = prior(3,1) + (prior(3,2) - prior(3,1)).*rand; % generate one random number within the given limits.
    % sigma_N (Variance noise floor)
    param_sigma_N(i) = prior(4,1) + (prior(4,2) - prior(4,1)).*rand; % generate one random number within the given limits.
    
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
        param_T(i);...
        param_G0(i);...
        param_lambda(i);...
        param_sigma_N(i)];
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

% Calculate first weights and covariance
weights_T = ones(1,nbr_extract)./nbr_extract;
weights_G0 = ones(1,nbr_extract)./nbr_extract;
weights_lambda = ones(1,nbr_extract)./nbr_extract;
weights_sigma_N = ones(1,nbr_extract)./nbr_extract;
var_T = var(params_T(1,:));
var_G0 = var(params_G0(1,:));
var_lambda = var(params_lambda(1,:));
var_sigma_N = var(params_sigma_N(1,:));

%% Iterations
p_theta = zeros(sumstat_iter,1);
for a = 2:iterations
    out = zeros(6,sumstat_iter);
    d = zeros(sumstat_iter,1);
    
    % Choose theta from accepted parameters of last iteration with propbability based on wieghts
    index = randsample((1:nbr_extract),sumstat_iter,true,weights);
    theta_prop = [params_T(a-1,index)' params_G0(a-1,index)' params_lambda(a-1,index)' params_sigma_N(a-1,index)'];
    % Select weights corresponding to each selected theta
    weight = weights(index);    
    parfor i = 1:sumstat_iter
        % Perturb theta 
        theta_curr = normrnd(theta_prop(i,:),covariance);
        while(check_params(theta_curr,prior)==2)
            theta_curr = mvnrnd(theta_prop(i,:),covariance);
        end
        % calculate probability that theta was generated
        p_theta(i) = weight(i)*ones(1,4)*mvnpdf(theta_curr./diag(covariance),theta_prop(i,:)./diag(covariance),diag([1 1 1 1]));
        
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
            param_T(i);...
            param_G0(i);...
            param_lambda(i);...
            param_sigma_N(i);...
            p_theta(i)];
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
    prop_theta = out(6,1:nbr_extract);
    
   % Update the prior for the next iteration
   [f_T,xi_T] = ksdensity(params_T(a,:));
   [f_G0,xi_G0] = ksdensity(params_G0(a,:));
   [f_lambda,xi_lambda] = ksdensity(params_lambda(a,:));
   [f_sigma_N,xi_sigma_N] = ksdensity(params_sigma_N(a,:));

    for j = 1:nbr_extract
        weights(j) = f_T/sum(prop_theta);
    end
    covariance = 2*diag(diag(cov(out(2:5,1:nbr_extract)')));
    disp(a);
end
disp('ABC algorithm finished... ')
toc
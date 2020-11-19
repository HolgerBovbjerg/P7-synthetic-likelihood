%% ABC implementation 2 - ABC REJECTION ALGORITHM:
clear
%% --- Global turin model simulation parameters ---------------------------------------------------
N  = 300;    % Number of different turin simulations.
Ns = 801;   % Number of time entries for each turin simulation. 
Bw = 4e9;   % Bandwidth (4Ghz).

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
load('Prior_data_large_prior_min_max_values.mat')

%% --- Generate "observed data" -----------------------------------------------------
load("Theta_true_values.mat") 
load('covariance_large_prior.mat')
% [covariance, ~] = find_cov_prior(prior);
% covariance = covariance/1000;
%%
% load('S_obs.mat')
load('S_obs_meas_data.mat')

mu_S_obs = mean(S_obs);     % Mean of the summary statistics 
Sigma_S_obs = cov(S_obs);     % Covariance of summary statistics
 
%% --- ABC rejection algorithm ---------------------------------------------------------------------
% Number of summary statistics sets to generate  
k = 10e3;    % Number of MCMC steps
thetas = zeros(4,k);
summ_stats = zeros(9,k);
dists = zeros(k,1);
accept = 0;
tol = 100e3;

disp('ABC algorithm computing, please wait... ')
tic

% T (Reverberation time):
param_T = prior(1,1) + (prior(1,2) - prior(1,1)).*rand; % generate one random number
% G0 (Reverberation gain)  
param_G0 = prior(2,1) + (prior(2,2) - prior(2,1)).*rand; % generate one random number within the given limits.
% lambda ()  
param_lambda = prior(3,1) + (prior(3,2) - prior(3,1)).*rand; % generate one random number within the given limits.
% sigma_N (Variance noise floor)
param_sigma_N = prior(4,1) + (prior(4,2) - prior(4,1)).*rand; % generate one random number within the given limits.

theta_curr = [param_T param_G0 param_lambda param_sigma_N];

thetas(:,1) = theta_curr';
%% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
[Pv, t] = sim_turin_matrix_gpu_w_delay(N, Bw, Ns, theta_curr,6e-9);

s_curr = create_statistics(Pv, t);
summ_stats(:,1) = s_curr;
%% STEP 3: calculate the difference between observed and simulated summary statistics 
% Mahalanobis distance see formular in document.
d_curr = (s_curr - mu_S_obs)/Sigma_S_obs * (s_curr - mu_S_obs)';

dists(1) = d_curr;

loglike_curr = -0.5/tol*d_curr^2; %Gaussian kernel weighting function

%% STEP 4: Run MCMC sampling
for j = 1:k
    cla reset
    theta_prop = mvnrnd(theta_curr,covariance);
    while(check_params(theta_prop,prior)==2)
        theta_prop = mvnrnd(theta_curr,covariance);
    end
    %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
    [Pv, t] = sim_turin_matrix_gpu_w_delay(N, Bw, Ns, theta_prop,6e-9);
    s_prop = create_statistics(Pv, t);
    %% STEP 3: calculate the difference between observed and simulated summary statistics 
    % Mahalanobis distance see formular in document.
    d_prop = (s_prop - mu_S_obs)/Sigma_S_obs * (s_prop - mu_S_obs)';
    
    loglike_prop = -0.5/tol*d_prop^2; %Gaussian kernel weighting function
    
        % MCMC ABC accept/reject step
    if (rand<exp(loglike_prop - loglike_curr))
        theta_curr = theta_prop;
        d_curr = d_prop;
        s_curr = s_prop;
        loglike_curr = loglike_prop;
        accept = accept+1;
    end
    thetas(:,j) = theta_curr;
    summ_stats(:,j) = s_curr';
    dists(j) = d_curr;
    
    real_time_plots(theta_true,thetas,j-1,accept,k,prior,[1 1 1 1],loglike_prop,1);
end 


disp('ABC algorithm finished... ')
toc
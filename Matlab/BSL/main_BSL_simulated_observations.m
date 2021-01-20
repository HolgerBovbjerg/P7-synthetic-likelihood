%%
clear all
load('Prior_simulated_observations_min_max_values.mat')

N = 300; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
%% Find first proposed Theta and covariance for proposal distribution
[~, theta_curr] = find_cov_prior(prior);
theta_start = theta_curr;
load('covariance_small_prior.mat')
covariance = covariance/17;

%% "Observed data for testing"
% load('S_obs_9_stats.mat')
[~, theta_true] = find_cov_prior(prior);
[Pv, t] = sim_turin_matrix_gpu_w_delay(2000, B, Ns, theta_true, 1e-8);
s_obs = create_statistics(Pv, t);
%%
k = 3500;    % Number of MCMC steps
L = 500;     % Number of statistics vectors used per likelihood.

accept = 0;
s_sim = zeros(L,9);
thetas = zeros(4,k);
thetas(:,1) = theta_curr';

tic
for i = 1:L
    [Pv, t] = sim_turin_matrix_gpu_w_delay(N, B, Ns, theta_curr, 1e-8);
    s_sim(i,:) = create_statistics(Pv, t);
end
toc

loglikelihood = synth_loglikelihood(s_obs,s_sim);

%% MCMC
for j = 2:k
    j
    cla reset
    % Draw a new proposed theta from multivariate normal distribution
    theta_prop = mvnrnd(theta_curr,covariance);
    i = 0;
    % Check to see if the proposed theta is within the prior range
    % if not keep drawing  
    while(check_params(theta_prop,prior)==2) 
        theta_prop = mvnrnd(theta_curr,covariance);
        i = i + 1
    end
    % Create statistics based on the proposed theta 
    parfor i = 1:L
        [Pv, t] = sim_turin_matrix_gpu_w_delay(N, B, Ns, theta_prop, 1e-8);
        s_sim(i,:) = create_statistics(Pv, t);
    end
    % Calculate new log-likelihood from statistics
    loglikelihoodnew = (synth_loglikelihood(s_obs,s_sim));
    
    % Compare new and old likelihood and accept based on Hasting-ratio
    if exp(loglikelihoodnew-loglikelihood) > rand
        loglikelihood = loglikelihoodnew;
        theta_curr = theta_prop;
        accept = accept+1;
    end
    % Keep the most recent theta
    thetas(:,j) = theta_curr';
    
    real_time_plots(theta_true,thetas,j-1,accept,k,prior,theta_prop,loglikelihoodnew,loglikelihood);
end

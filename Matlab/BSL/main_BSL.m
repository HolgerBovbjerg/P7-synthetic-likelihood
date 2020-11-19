%%
clear all
load('Prior_data_large_prior_min_max_values.mat')
load('Theta_true_values.mat')

N = 300; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
%% Find first proposed Theta and covariance for proposal distribution
[covariance, theta_curr] = find_cov_prior(prior);
% theta_curr = theta_true;
theta_curr(2) = db2pow(-88);
theta_curr(3) = 1.89e9;
theta_start = theta_curr;
load('covariance_small_prior.mat')
% covariance = covariance/14;

%% "Observed data for testing"
load('S_obs_9_stats.mat')
% load('observed_data_statistics.mat')
%%
k = 2500;    % Number of MCMC steps
L = 400;     % Numberof statistics vectors used per likelihood.

accept = 0;
s_sim = zeros(L,9);
thetas = zeros(4,k);
thetas(:,1) = theta_curr';
tic
parfor i = 1:L
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr);
    s_sim(i,:) = create_statistics(Pv, t);
    i
end
toc
loglikelihood = synth_loglikelihood(s_obs,s_sim);

%% MCMC

for j = 108:k
    j
    
    if accept/j < 0.25 && mod(j,100)==0
        covariance = covariance/1.1;
        loglikelihood = loglikelihoodnew
        if mod(j,150)==0
            L = L + 10
        end
    end
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
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop);
        s_sim(i,:) = create_statistics(Pv, t);
    end
    % Calculate new log-likelihood from statistics
    loglikelihoodnew = (synth_loglikelihood(s_obs,s_sim));
    
    % Compare new and old likelihood and accept based on Hasting-ratio
    if exp(loglikelihoodnew-loglikelihood) >  rand
        loglikelihood = loglikelihoodnew;
        theta_curr = theta_prop;
        accept = accept+1;
    end
    % Keep the most recent theta
    thetas(:,j) = theta_curr';
    
    real_time_plots(theta_true,thetas,j-1,accept,k,prior,theta_prop,loglikelihoodnew,loglikelihood);
end

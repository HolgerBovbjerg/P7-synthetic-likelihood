%%
clear all
load('Prior_data_very_small_prior_min_max_values.mat')
load('Theta_true_values.mat')

N = 50; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
%%
[covariance, theta_curr] = find_cov_prior(prior);
covariance = covariance/100;
theta_start = theta_curr;
%%
% "Observed data for testing"
load('S_obs_true.mat')

% [Pv, t] = sim_turin_matrix_gpu(10000, B, Ns, theta_true);
% s_obs = create_statistics(Pv, t);
%%
L = 10;     % Numberof statistics vectors used per likelihood
k = 500;    % Number of MCMC steps
accept = 0;
s_sim = zeros(L,4);
thetas = zeros(4,k);
thetas(:,1) = theta_curr';

parfor i = 1:L
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr);
    s_sim(i,:) = create_statistics(Pv, t);
end

loglikelihood = synth_loglikelihood(s_obs,s_sim);

for j = 2:k
    cla reset
    theta_prop = mvnrnd(theta_curr,covariance);
    i = 0;
    while(check_params(theta_prop,prior)==2)
        theta_prop = mvnrnd(theta_curr,covariance);
        i = i + 1
    end
    
    parfor i = 1:L
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop);
        s_sim(i,:) = create_statistics(Pv, t);
    end
    loglikelihoodnew = (synth_loglikelihood(s_obs,s_sim));
    
    if exp(loglikelihoodnew-loglikelihood) >  rand
        loglikelihood = loglikelihoodnew;
        theta_curr = theta_prop;
        accept = accept+1;
    end
    thetas(:,j) = theta_curr';
    
    real_time_plots(theta_true,thetas,j-1,accept,k,prior);
end

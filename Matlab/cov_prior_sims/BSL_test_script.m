%%
clear all
load('Prior_data_large_prior_min_max_values.mat')
load('Theta_true_values.mat')

N = 10; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
%%
[covariance, theta_curr] = find_cov_prior(prior);
theta_start = theta_curr;
load("covariance_large_prior.mat");

scale = 1/10;
covariance = covariance*scale;
%%
% "Observed data for testing"
% load('S_obs_true.mat')

[Pv, t] = sim_turin_matrix_gpu(1000, B, Ns, theta_true);
s_obs = create_statistics_new(Pv, t);
%%
k = 2000;    % Number of MCMC steps
L = 200;     % Numberof statistics vectors used per likelihood

accept = 0;
s_sim = zeros(L,18);
thetas = zeros(4,k);
thetas(:,1) = theta_curr';

parfor i = 1:L
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr);
    s_sim(i,:) = create_statistics_new(Pv, t);
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
        s_sim(i,:) = create_statistics_new(Pv, t);
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

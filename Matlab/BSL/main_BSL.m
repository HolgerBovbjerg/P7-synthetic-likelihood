%%
clear all
load('Prior_data_large_prior_min_max_values.mat')
load('Theta_true_values.mat')

N = 300; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
%% Find first proposed Theta and covariance for proposal distribution
[covariance, theta_curr] = find_cov_prior(prior);

theta_start = theta_curr;
load('covariance_medium_prior.mat')

%% "Observed data for testing"
load('S_obs_9_stats.mat')
%%
k = 20000;    % Number of MCMC steps
L = 500;     % Numberof statistics vectors used per likelihood

accept = 0;
s_sim = zeros(L,9);
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
    
    real_time_plots(theta_true,thetas,j-1,accept,k,prior,theta_prop,loglikelihoodnew,loglikelihood);
end

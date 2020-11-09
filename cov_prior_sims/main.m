%%
clear all
N = 200; % Number of Turin simulations

T_true = 7.8e-9;
G0_true = db2pow(-83.9); % Reverberation gain converted from dB to power - 4.0738e-9
lambda_true = 10e7;
sigma_N_true = sqrt(0.28e-9); % noise variance

theta_true = [T_true G0_true lambda_true sigma_N_true];
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
%%
[covariance theta_curr] = find_cov_prior(G0_true);
covariance = covariance/100
theta_start = theta_curr;
%%
% -------------------------------------------- %
% "Observed data for testing"

[Pv, t] = sim_turin_matrix_gpu(N*10, B, Ns, theta_true);
s_obs = create_statistics(Pv, t);

L = 10; % Numberof statistics vectors used per likelihood

s_sim = zeros(L,8); 
for i = 1:L
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr);
    s_sim(i,:) = create_statistics(Pv, t);
end
  
loglikelihood = synth_loglikelihood(s_obs,s_sim)

accept = 0;
 
k = 500;    % Number of MCMC steps
thetas = zeros(4,k);
thetas(:,1) = theta_curr';

tic
for j = 2:k
    cla reset
    j
    theta_prop = mvnrnd(theta_curr,covariance);
    while(any(theta_prop < 0))
    theta_prop = mvnrnd(theta_curr,covariance);
    end
covs_T(j) = covariance(1,1);
covs_G0(j) = covariance(2,2);
covs_lambda(j) = covariance(3,3);
covs_sigma(j) = covariance(4,4);
    parfor i = 1:L
      [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop);
      s_sim(i,:) = create_statistics(Pv, t);
    end
    loglikelihoodnew = (synth_loglikelihood(s_obs,s_sim));

    likeli = exp(loglikelihoodnew-loglikelihood);
    if likeli >  rand 
    loglikelihood = loglikelihoodnew;
    theta_curr = theta_prop;
    accept = accept+1
    end
    thetas(:,j) = theta_curr';
    
real_time_plots(thetas,j,accept,k);
end   
toc

  %%
  figure(1)

%%
clear all
N = 200; % Number of Turin simulations

T_true = 7.8e-9;
G0_true = db2pow(-83.9); % Reverberation gain converted from dB to power - 4.0738e-9
lambda_true = 10e8;
sigma_N_true = sqrt(0.28e-9); % noise variance

theta_true = [T_true G0_true lambda_true sigma_N_true];
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
%%
[covariance theta_curr] = find_cov_prior;
theta_start = theta_curr;
the
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
h = 10000;  % Number of MCMC steps
k = 200;
thetas = zeros(4,h+k);
likelihoods = zeros(1,N);
thetas(:,1) = theta_curr';
covs_T = zeros(1,h+k);
covs_G0 = zeros(1,h+k);
covs_lambda = zeros(1,h+k);
covs_sigma = zeros(1,h+k);
tic
for j = 2:k
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
end   
toc

%%
k=j
for j = k+1:h+k
    j
    covariance = adaptive_cov(thetas(:,1:j-1),covariance,j-1);
covs_T(j) = covariance(1,1);
covs_G0(j) = covariance(2,2);
covs_lambda(j) = covariance(3,3);
covs_sigma(j) = covariance(4,4);
    theta_prop = mvnrnd(theta_curr,covariance);
    while(any(theta_prop < 0))
    theta_prop = mvnrnd(theta_curr,covariance);
    end
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
end
  %%
tiledlayout(4,1)
nexttile
plot(thetas(1,1:j),'o')
yline(T_true)
title("T")

nexttile
plot(thetas(2,1:j),'o')
yline(G0_true)
title("G0")

nexttile
plot(thetas(3,1:j),'o')
yline(lambda_true)
title("\lambda")

nexttile
plot(thetas(4,1:j),'o')
yline(sigma_N_true)
title("\sigma")
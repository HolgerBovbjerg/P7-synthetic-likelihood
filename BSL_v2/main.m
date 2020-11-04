clear all


N = 200; % Number of Turin simulations

T_true = 7.8e-9;
lambda_true = 10e8;

G0_true = db2pow(-83.9); % Reverberation gain converted from dB to power - 4.0738e-9
sigma_N_true = sqrt(0.28e-9); % noise variance
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

[covariance theta_curr] = find_cov_prior;

% -------------------------------------------- %
% "Observed data for testing"

[Pv, t] = sim_turin_matrix_gpu(N*10, B, Ns, T_true, G0_true, lambda_true, sigma_N_true);
s_obs = create_statistics(Pv, t);

L = 10; % Numberof statistics vectors used per likelihood

s_sim = zeros(L,8); 
for i = 1:L
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr(1), theta_curr(2), theta_curr(3), theta_curr(4));
    s_sim(i,:) = create_statistics(Pv, t);
end
  
loglikelihood = synth_loglikelihood(s_obs,s_sim)

accept = 0;
h = 120;
k = 100; % Number of MCMC steps
thetas= zeros(4,h);
likelihoods = zeros(1,N);
thetas(:,1) = theta_curr';
  for j = 2:k
      j
      L = 10; % Numberof statistics vectors used per likelihood
      while(any(theta_prop < 0))
        theta_prop = mvnrnd(theta_curr,covariance);
      end

      for i = 1:L
          [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop(1), theta_prop(2), theta_prop(3), theta_prop(4));
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

  for j = k+1:h
      j
      covariance = adaptive_cov(thetas(:,1:j-1),covariance,j-1);
      while(any(theta_prop < 0))
        theta_prop = mvnrnd(theta_curr,covariance);
      end
      for i = 1:L
          [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop(1), theta_prop(2), theta_prop(3), theta_prop(4));
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
clear all

N = 200; % Number of Turin simulations
likelihoods = zeros(1,N);
T = 7.8e-9;
lambda = 10e8;
B = 4e9; % Bandwidth of signal: 4 GHz
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power - 4.0738e-9

sigma_N = sqrt(0.28e-9); % noise variance
Ns = 801; % Number of sample points per Turin simulation

theta_para_cov = find_theta_cov;
% -------------------------------------------- %
% "Observed data for testing"

[Pv, t] = sim_turin_matrix_gpu(N*10, B, Ns, T, G0, lambda, sigma_N);

s_obs = create_statistics(Pv, t);

% -------------------------------------------- %

% Initial guess - first value of the chain 
T =5.055e-9;
tic
lambda = 0.599878e+09;
G0 = 6.056e-09; % Reverberation gain converted from dB to power
sigma_N = 1.073155e-05 ; % noise variance
L = 10; % Numberof statistics vectors used per likelihood

theta_curr = [T G0 lambda sigma_N];
s_sim = zeros(L,8); 
for i = 1:L
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr(1), theta_curr(2), theta_curr(3), theta_curr(4));
    s_sim(i,:) = create_statistics(Pv, t);
end
  
 theta_mean = mean(s_sim);
 theta_cov = cov(s_sim)
 %theta_cov2 = (1/(length(s_obs)-1))*sum(((s_obs-theta_mean)*(s_obs-theta_mean)'))
 loglikelihood = (synth_loglikelihood(s_obs,s_sim))
 % -------------------------------------------- %
  %     % MCMC part
  accept = 0;
  k = 3000; % Number of MCMC steps
  thetas= zeros(4,k);
  
  for j = 1:k
      j
      L = 10; % Numberof statistics vectors used per likelihood

         % theta_prop = mvnrnd(theta_curr,theta_para_cov);

      

          theta_prop(1) =normrnd(theta_curr(1),sqrt(theta_para_cov(1,1)));
          theta_prop(2) =normrnd(theta_curr(2),sqrt(theta_para_cov(2,2)));
          theta_prop(3) =normrnd(theta_curr(3),sqrt(theta_para_cov(3,3)));
          theta_prop(4) =normrnd(theta_curr(4),sqrt(theta_para_cov(4,4)));

      for i = 1:L
          [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop(1), theta_prop(2), theta_prop(3), theta_prop(4));
          s_sim(i,:) = create_statistics(Pv, t);
      end
      loglikelihoodnew = (synth_loglikelihood(s_obs,s_sim));

      if exp(loglikelihoodnew-loglikelihood) > rand
%        if r > rand
        loglikelihood = loglikelihoodnew;
          theta_curr = theta_prop;
          accept = accept+1
      end
      thetas(:,j) = theta_curr'; 
  end   
  toc
  %%
    figure(2)
    tiledlayout(4,1)
    nexttile
    
    plot(1:k,thetas(1,:),"o")
    yline(7.8e-9)
   % ylim([6e-9 9e-9])
    title("T")
    
    nexttile
    
    plot(1:k,thetas(2,:),"o")
    yline(db2pow(-83.9))
   % ylim([3e-9 6e-9])
    title("G0")
    
    nexttile
    
    plot(1:k,thetas(3,:),"o")
    yline(10e8)
   % ylim([0.5e+9 1.5e+9])
    title("\lambda")
    
    nexttile
    
    plot(1:k,thetas(4,:),"o")
    yline(sqrt(0.28e-9))
    %ylim([1.0e-5 2.5e-5])
    title("\sigma")
    
%   %%
%     figure(3)
%     tiledlayout(4,1)
%     nexttile
%     
%     plot(1:1000,theta_prop(1,:),"o")
%     yline(7.8e-9)
%     ylim([6e-9 9e-9])
%     title("T")
%     
%     nexttile
%     
%     plot(1:1000,theta_prop(2,:),"o")
%     yline(db2pow(-83.9))
%     ylim([3.5e-9 5.2e-9])
%     title("G0")
%     
%     nexttile
%     
%     plot(1:1000,theta_prop(3,:),"o")
%     yline(10e8)
%     ylim([0.85e+9 1.4e+9])
%     title("\lambda")
%     
%     nexttile
%     
%     plot(1:1000,theta_prop(4,:),"o")
%     yline(sqrt(0.28e-9))
%     ylim([1.4e-5 2.1e-5])
%     title("\sigma")
%     
%       %%
%     figure(4)
%     tiledlayout(4,1)
%     nexttile
%     
%     plot(1:1000,theta_prop2(1,:),"o")
%     yline(7.8e-9)
%     ylim([6e-9 9e-9])
%     title("T")
%     
%     nexttile
%     
%     plot(1:1000,theta_prop2(2,:),"o")
%     yline(db2pow(-83.9))
%     ylim([3.5e-9 5.2e-9])
%     title("G0")
%     
%     nexttile
%     
%     plot(1:1000,theta_prop2(3,:),"o")
%     yline(10e8)
%     ylim([0.85e+9 1.4e+9])
%     title("\lambda")
%     
%     nexttile
%     
%     plot(1:1000,theta_prop2(4,:),"o")
%     yline(sqrt(0.28e-9))
%     ylim([1.4e-5 2.1e-5])
%     title("\sigma")
%     



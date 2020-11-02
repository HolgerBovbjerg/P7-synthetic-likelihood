clear all

N = 200; % Number of Turin simulations

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
parfor i = 1:L
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr(1), theta_curr(2), theta_curr(3), theta_curr(4));
    s_sim(i,:) = create_statistics(Pv, t);
end
  
 theta_mean = mean(s_sim);
 theta_cov = cov(s_sim)
 %theta_cov2 = (1/(length(s_obs)-1))*sum(((s_obs-theta_mean)*(s_obs-theta_mean)'))
 loglikelihood = synth_loglikelihood(s_obs,s_sim)
 % -------------------------------------------- %
  %     % MCMC part
  accept = 0;
  k = 200; % Number of MCMC steps
  thetas= zeros(4,k);
likelihoods = zeros(1,N);
  thetas(:,1) = theta_curr'
  for j = 2:k
      j
      L = 10; % Numberof statistics vectors used per likelihood

         % theta_prop = mvnrnd(theta_curr,theta_para_cov);

      

          theta_prop(1) =abs(normrnd(theta_curr(1),sqrt(theta_para_cov(1,1))));
          theta_prop(2) =abs(normrnd(theta_curr(2),sqrt(theta_para_cov(2,2))));
          theta_prop(3) =abs(normrnd(theta_curr(3),sqrt(theta_para_cov(3,3))));
          theta_prop(4) =abs(normrnd(theta_curr(4),sqrt(theta_para_cov(4,4))));

      parfor i = 1:L
          [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop(1), theta_prop(2), theta_prop(3), theta_prop(4));
          s_sim(i,:) = create_statistics(Pv, t);
      end
      loglikelihoodnew = (synth_loglikelihood(s_obs,s_sim));
    
    likeli = exp(loglikelihoodnew-loglikelihood);
      if likeli >  rand 
%        if r > rand
        loglikelihood = loglikelihoodnew;
          theta_curr = theta_prop;
          accept = accept+1
      end
      thetas(:,j) = theta_curr';
      likelihoods(j) = loglikelihoodnew;
  end   
  toc
  %%
  iters = 2;
  values = [0 0 0 0]';
for q = 1:iters
    values = [values thetas];
    pd1 = fitdist(thetas(1,:)','Normal');
    pd2 = fitdist(thetas(2,:)','Normal');
    pd3 = fitdist(thetas(3,:)','Normal');
    pd4 = fitdist(thetas(4,:)','Normal');
    sigma1 = pd1.sigma;
    sigma2 = pd2.sigma;
    sigma3 = pd3.sigma;
    sigma4 = pd4.sigma;
    theta_curr(1) = pd1.mu;
    theta_curr(2) = pd2.mu;
    theta_curr(3) = pd3.mu;
    theta_curr(4) = pd4.mu;
    
    
    tic
    N = 200; % Number of Turin simulations
    B = 4e9; % Bandwidth of signal: 4 GHz
    Ns = 801;
    
    L = 10;
    s_sim = zeros(L,8);
    parfor i = 1:L
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr(1), theta_curr(2), theta_curr(3), theta_curr(4));
        s_sim(i,:) = create_statistics(Pv, t);
    end
    loglikelihood = synth_loglikelihood(s_obs,s_sim)/100
    
    %theta_para_cov = find_theta_cov;
    
    
    clear thetas
    accept = 0;
    k = 200; % Number of MCMC steps
    thetas= zeros(4,k);
    likelihoods = zeros(1,N);
    thetas(:,1) = theta_curr'
    for j = 2:k
        j
        L = 10; % Numberof statistics vectors used per likelihood
        
        % theta_prop = mvnrnd(theta_curr,theta_para_cov);
        
        
        
        theta_prop(1) =abs(normrnd(theta_curr(1),sigma1));
        theta_prop(2) =abs(normrnd(theta_curr(2),sigma2));
        theta_prop(3) =abs(normrnd(theta_curr(3),sigma3));
        theta_prop(4) =abs(normrnd(theta_curr(4),sigma4));
        
        parfor i = 1:L
            [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop(1), theta_prop(2), theta_prop(3), theta_prop(4));
            s_sim(i,:) = create_statistics(Pv, t);
        end
        loglikelihoodnew = (synth_loglikelihood(s_obs,s_sim));
        
        likeli = exp(loglikelihoodnew-loglikelihood);
        if likeli >  rand
            %        if r > rand
            loglikelihood = loglikelihoodnew;
            theta_curr = theta_prop;
            accept = accept+1
        end
        thetas(:,j) = theta_curr';
        likelihoods(j) = likeli;
    end
end
   toc
   
   %%
   
      figure(5)
    tiledlayout(4,1)
    nexttile
    
    plot(1:length(values),values(1,:),"o")
    yline(7.8e-9)
   % ylim([6e-9 9e-9])
    title("T")
    
    nexttile
    
    plot(1:length(values),values(2,:),"o")
    yline(db2pow(-83.9))
   % ylim([3e-9 6e-9])
    title("G0")
    
    nexttile
    
    plot(1:length(values),values(3,:),"o")
    yline(10e8)
   % ylim([0.5e+9 1.5e+9])
    title("\lambda")
    
    nexttile
    
    plot(1:length(values),values(4,:),"o")
    yline(sqrt(0.28e-9))
    %ylim([1.0e-5 2.5e-5])
    title("\sigma")
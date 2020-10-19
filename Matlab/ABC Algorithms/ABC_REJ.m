
%% ABC implementation 1 most basic - ABC REJECTION ALGORITHM:
% Population Monte Carlo simulation 
% based om the paper: 
% Aproximate Bayesian Computation usink Markov Chain Monte Carlo simulation: DREAM(ABC) 

% Once summary statistics have been defined we are left with finding those
% values of the parameter for which <= epsilon.


%% ------ Generate "observed data" used as Y_obs -----------------------------

param_T       = 7.8e-12; % -> 
param_G0      = -19.24;  % -> 
param_lambda  = 10e-12;  % -> 
sigma_N       = 2.8e-10; % ->

N  = 50;        % Number of turin simulations (Does not affect the size of OUTPUT P_Y)
Ns = 650;       % Number of Hk's (transfer functions)
Bw = 5000000;   % Bandwidth (5Mhz)
[P_Y_obs, t_obs] = sim_turin_matrix(N, Bw, Ns, param_T, param_G0, param_lambda, sigma_N);

nth_moment = 3; 
% Do summary statistics on the observed dataset
[X, S_Y_obs] = sumStat(t_obs,P_Y_obs, nth_moment); 

% -------- END "observed data"  -------------------------------------------
%%

iter = 1; % Set the number of simulations to peform
epsilon = 1e-40; % Acceptable error 
dist_func = epsilon + 1;
accepted_parameters = zeros(1,iter);

dist_func_vec = zeros(1000,1);
tic
%% ABC Rejection algorith 
for i = 1:iter
%  OBS -> Need to do several for loops in order to test all three parameters
%  T, G0 and lambda
    while dist_func > epsilon 
       
        %% STEP 1: 
        % Sample parameter from the prior (for this we need a prior distribution?) 
        % Need to be pulled from a distribution
       
         % T (Reverberation time) sample from prior:
         T_a = 7.8e-11; 
         T_b = 7.8e-7;
         param_T = T_a + (T_b-T_a).*rand(1,1); % generate one random number
         
         % G0 (Reverberation gain) sample from prior:
         G0_a = -140; 
         G0_b = -70;
         param_G0 = db2pow(G0_a + (G0_b-G0_a).*rand(1,1)); % generate one random number
     
         % lambda () sample from prior:
         lambda_a = 10e8; 
         lambda_b = 10e10;
         param_lambda = lambda_a + (lambda_b-lambda_a).*rand(1,1); % generate one random number
         
         % sigma_N (Variance noise floor) sample from prior:
         sigmaN_a = sqrt(0.28e-10); 
         sigmaN_b = sqrt(0.28e-8);
         param_sigma_N = sigmaN_a + (sigmaN_b-sigmaN_a).*rand(1,1); % generate one random number
        
         %% STEP 2:
         % Simulate data using Turing model based on sample parameter from step 1 
         N = 50;        % Number of turin simulations (Does not affect the size of OUTPUT P_Y)
         Ns = 650;      % Number of Hk's (transfer functions)
         Bw = 4e9;      % Bandwidth (4Ghz)
        
         % P_Y is a row vector of Ns length with the mean of N turin
         % simulations
         [P_Y, t] = sim_turin_matrix(N, Bw, Ns, param_T, param_G0, param_lambda, param_sigma_N);
        
        
         %% STEP 3:
         %  Do summary statistics on the simulated dataset: 
         nth_moment = 3; 
         [S_X, S_X_mean] = sumStat(t, P_Y, nth_moment); 
      
         %% STEP 4:
         % Calculate the distance function (difference between simulated and observed summary statistics)  
         % sum_stat_difference = ABC_formula     
         dist_func = abs(S_Y_obs) - S_X_mean;
         dist_func_vec(i) = dist_func;
    end  
end 
toc

xx = 1:1000;
plot(xx,dist_func_vec)

% 
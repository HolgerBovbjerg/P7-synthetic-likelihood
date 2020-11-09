%% ABC implementation 2 - ABC REJECTION ALGORITHM:
clear
%% --- Global turin model simulation parameters ---------------------------------------------------
N  = 50;    % Number of different turin simulations.
Ns = 801;   % Number of time entries for each turin simulation. 
Bw = 4e9;   % Bandwidth (4Ghz).

%% --- Generate "observed data" -----------------------------------------------------
param_T_obs       = 7.8e-9;  % Reverberation time
param_G0_obs      = db2pow(-83.9);    % linear gain 
param_lambda_obs  = 10e9;       % arrival rate (1/s)    
param_sigma_N_obs = sqrt(0.28e-9);   % Noise std
M = 200; %number of summary statisctics realisations

% S_obs = create_statistics(M, N, param_T_obs , param_G0_obs,...
% param_lambda_obs, param_sigma_N_obs, Bw, Ns); % create new summary statistic of
% observed data
load S_obs; % Load already generated data

%%
mu_S_obs = mean(S_obs);     % Mean of the mean and varaince of log(moments)?
Sigma_S_obs = cov(S_obs);   

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
% a = min , b = max
 T_a = 1e-9; 
 T_b = 15e-9;  
 G0_a = db2pow(pow2db(param_G0_obs) - 10);    % Power gain (not in dB)
 G0_b = db2pow(pow2db(param_G0_obs) + 10);     % Power gain (not in dB)
 lambda_a = 1e8;
 lambda_b = 20e9;
 sigmaN_a = sqrt(0.28e-10); 
 sigmaN_b = sqrt(0.28e-8);
 
%% --- ABC rejection algorithm ---------------------------------------------------------------------
% Set total iterations
iterations = 3;

% Number of summary statistics sets to generate  
sumstat_iter = 2000;

% Extract this amount of parameter entries from each generated summary
% statistic
nbr_extract = 100;

% Preallocation of vectors: 
out = zeros(5,sumstat_iter);
meanVar_params = zeros(8,iterations);

params_T = zeros(iterations,nbr_extract);
params_G0 = zeros(iterations,nbr_extract);
params_lambda = zeros(iterations,nbr_extract);
params_sigma_N = zeros(iterations,nbr_extract);

disp('ABC algorithm computing, please wait... ')
tic
for a = 1:iterations 
    out = zeros(5,sumstat_iter);
    d = zeros(sumstat_iter,1);
    parfor i = 1:sumstat_iter
        %% STEP 1: Sample parameter from predefined prior distribution (uniform):      
        % T (Reverberation time):
        param_T = T_a + (T_b-T_a)*rand; % generate one random number
        % G0 (Reverberation gain)  
        param_G0 = (G0_a + (G0_b-G0_a)*rand); % generate one random number within the given limits.
        % lambda ()  
        param_lambda = lambda_a + (lambda_b-lambda_a)*rand; % generate one random number within the given limits.
        % sigma_N (Variance noise floor)
        param_sigma_N = sigmaN_a + (sigmaN_b-sigmaN_a)*rand; % generate one random number within the given limits.
       
        %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
        cd ../        % change folder for statistics function
        cd statistics
        S_simulated = create_statistics(1, N, Bw, Ns, 'matrix', param_T, param_G0, param_lambda, param_sigma_N);
        %% STEP 3: calculate the difference between observed and simulated summary statistics 
        % Mahalanobis distance see formular in document.
        d(i) = (S_simulated - mu_S_obs)/Sigma_S_obs * (S_simulated - mu_S_obs)';
        
        % Row 1 of the out vector contains the distance 
        % the rest of the rows contains the corresponding parameters 
        % used for generating that specific distance.    
        out(:,i) =  [d(i);...
                    param_T;...
                    param_G0;...
                    param_lambda;...
                    param_sigma_N];
        disp(i);
    end 
    % Sort the "out" matrix so that the lowest euclidean distance is at the
    % (1,1) matrix position and highest distance is at (max,1) 
    out = sortrows(out',1)';

    % Following vectors holds ALL extracted parameter values 
    % that was within euclidean distance
    params_T(a,:)       = out(2,1:nbr_extract);
    params_G0(a,:)      = out(3,1:nbr_extract);
    params_lambda(a,:)  = out(4,1:nbr_extract);
    params_sigma_N(a,:) = out(5,1:nbr_extract);

    
   % Update the min/max values for the parameters for the next iteration 
   T_a      = min(params_T(a,:));
   T_b      = max(params_T(a,:));
   G0_a     = min(params_G0(a,:));
   G0_b     = max(params_G0(a,:));
   lambda_a = min(params_lambda(a,:));
   lambda_b = max(params_lambda(a,:));
   sigmaN_a = min(params_sigma_N(a,:));
   sigmaN_b = max(params_sigma_N(a,:));
   disp(a);
end 
toc
   



 
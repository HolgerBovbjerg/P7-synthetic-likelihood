%% ABC implementation 2 - ABC REJECTION ALGORITHM:
% IMPORTANT NOTE: Need to be started with MATLAB in the folder ABC_algorithms to run withou errors
clear
%% --- Global turin model simulation parameters ---------------------------------------------------
N  = 50;    % Number of different turin simulations.
Ns = 801;   % Number of time entries for each turin simulation. 
Bw = 4e9;   % Bandwidth (4Ghz).

%% --- Generate "observed data" -----------------------------------------------------
T_obs           = 7.8e-9;           % Reverberation time
G0_obs          = db2pow(-83.9);    % linear gain 
lambda_obs      = 10e9;             % arrival rate (1/s)    
sigma_N_obs     = sqrt(0.28e-9);    % Noise std
M = 10;                             % Number of summary statisctics realisations

cd ../         % Change folder for accessing create statistics function
cd statistics
% Generate observed data summary statistics based on turin model and temporal moments
S_obs = create_statistics(M, N, Bw, Ns, 'matrix', T_obs , G0_obs, lambda_obs, sigma_N_obs);

% load S_obs;   % Uncomment to use

mu_S_obs = mean(S_obs);     % Mean of all the summary statistic vector   
Sigma_S_obs = cov(S_obs);   % covariance matrix of all the summary statustic vectors

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
% a = min, b = max
 T_a = 15e-9; 
 T_b = 1e-9;  
 G0_a = db2pow(-70);     % Power gain
 G0_b = db2pow(-90);     % Power gain
 lambda_a = 20e8;
 lambda_b = 20e9;
 sigmaN_a = sqrt(0.5e-9); 
 sigmaN_b = sqrt(0.1e-9);
 
%% --- ABC rejection algorithm ---------------------------------------------------------------------
% Set total iterations
iterations = 6;

% Number of summary statistics sets to generate  
sumstat_iter = 20;

% Extract this amount of parameter entries from each generated summary
% statistic
nbr_extract = 10;

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
    for i = 1:sumstat_iter
        %% STEP 1: Sample parameter from predefined prior distribution (uniform):   
        % Generate random numbers within the given limits.
        % T (Reverberation time):
        param_T = T_a + (T_b-T_a)*rand;
        % G0 (Reverberation gain)  
        param_G0 = (G0_a + (G0_b-G0_a)*rand);
        % lambda (Arrival rate)  
        param_lambda = lambda_a + (lambda_b-lambda_a)*rand; 
        % sigma_N (Variance noise floor)
        param_sigma_N = sigmaN_a + (sigmaN_b-sigmaN_a)*rand;
        
        cd ../        % Change folder for accessing create statistics function
        cd statistics
        %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
        S_simulated = create_statistics(1, N, Bw, Ns, 'matrix', param_T , param_G0, param_lambda, param_sigma_N);
        %% STEP 3: calculate the difference between observed and simulated summary statistics 
        % Mahalanobis distance see formular in worksheet.
        d(i) = (S_simulated - mu_S_obs)/Sigma_S_obs * (S_simulated - mu_S_obs)';
        %d(i) = mahal(S_simulated, S_obs);
        
        % Row 1 of the out vector contains the distance 
        % the rest of the rows contains the corresponding parameters 
        % used for generating that specific distance.    
        out(1,i) =  d(i);         
        out(2,i) =  param_T;
        out(3,i) =  param_G0;
        out(4,i) =  param_lambda;
        out(5,i) =  param_sigma_N;
        
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

% Plot distribution based after each iteration for G_0

t = tiledlayout(3,2, 'TileSpacing', 'none', 'Padding', 'compact');
title(t,'Parameter estimation using ABC rejection algorithm for Turin model');

nexttile
 % Plot the posterior distribution for the T parameter: 
 vect1 = params_G0(1,:);
 [f1,x1] = ksdensity(vect1);
 plot(x1,f1);
 xline(T_obs); % Original value
 title('T_1');

nexttile 
 % Plot the posterior distribution for the T parameter: 
 vect2 = params_G0(2,:);
 [f2,x2] = ksdensity(vect2);
 plot(x2,f2);
 xline(T_obs); % Original value
 title('T_2');
 
nexttile 
 % Plot the posterior distribution for the T parameter: 
 vect3 = params_G0(3,:);
 [f3,x3] = ksdensity(vect3);
 plot(x3,f3);
 xline(T_obs); % Original value
 title('T_3');

nexttile 
 % Plot the posterior distribution for the T parameter: 
 vect4 = params_G0(4,:);
 [f4,x4] = ksdensity(vect4);
 plot(x4,f4);
 xline(T_obs); % Original value
 title('T_4');
 
nexttile
 % Plot the posterior distribution for the T parameter: 
 vect5 = params_G0(5,:);
 [f5,x5] = ksdensity(vect5);
 plot(x5,f5);
 xline(T_obs); % Original value
 title('T_5');
 





 
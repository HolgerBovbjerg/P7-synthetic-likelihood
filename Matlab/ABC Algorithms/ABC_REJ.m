
%% ABC implementation 2 - ABC REJECTION ALGORITHM:

%% --- Global turin model simulation parameters ---------------------------------------------------
N  = 50;    % Number of different turin simulations.
Ns = 600;   % Number of time entries for each turin simulation. 
Bw = 4e9;   % Bandwidth (4Ghz).

%% --- Generate "observed data" -----------------------------------------------------
param_T       = 7.8e-9; 
param_G0      = 4.07e-9;    % dB = -83.9  
param_lambda  = 10e9;       % arrival rate (1/s)    
sigma_N       = 1.673e-4;   % equal to sqrt(28e-9)

[P_Y_observed, t_observed] = sim_turin_matrix(N, Bw, Ns, param_T, param_G0, param_lambda, sigma_N);

% Do summary statistics on the observed data 
S_observed = sumStatMeanMoment(t_observed, P_Y_observed); 
disp('Summary statistics of observed data generated...')

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
% a = min , b = max
 T_a = 7.8e-11; 
 T_b = 7.8e-7;  
 G0_a = 8.89e-11;    % Power gain (not in dB)
 G0_b = 4.07e-7;     % Power gain (not in dB)
 lambda_a = 1e8;
 lambda_b = 20e10;
 sigmaN_a = 1.673e-6; 
 sigmaN_b = 1.673e-2;
 
%% --- ABC rejection algorithm ---------------------------------------------------------------------
% Rejection ratio = 1 - (nbr_extract / sumstat_iter) 

% Set total iterations
iterations = 100;

% Number of summary statistics sets to generate  
sumstat_iter = 200;

% Extract this amount of parameter entries from each generated summary
% statistic
nbr_extract = 20;

% Preallocation of vectors: 
out = zeros(5,sumstat_iter);
meanVar_params = zeros(8,iterations);

index = 1; % used as index for vector holding ALL extracted parameter values 

disp('ABC algorithm computing, please wait... ')
tic
for a = 1:iterations 
    for i = 1:sumstat_iter
        %% STEP 1: Sample parameter from predefined prior distribution (uniform):
        
        % T (Reverberation time):
        param_T = T_a + (T_b-T_a).*rand(1,1); % generate one random number
        % G0 (Reverberation gain)  
        param_G0 = (G0_a + (G0_b-G0_a).*rand(1,1)); % generate one random number within the given limits.
        % lambda ()  
        param_lambda = lambda_a + (lambda_b-lambda_a).*rand(1,1); % generate one random number within the given limits.
        % sigma_N (Variance noise floor)
        param_sigma_N = sigmaN_a + (sigmaN_b-sigmaN_a).*rand(1,1); % generate one random number within the given limits.
        
        %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 
        [P_Y, t] = sim_turin_matrix(N, Bw, Ns, param_T, param_G0, param_lambda, param_sigma_N);
        
        %% STEP 3: Do summary statistics on the simulated dataset: 
        S_simulated = sumStatMeanMoment(t, P_Y);
         
        %% STEP 4: calculate the difference between observed and simulated summary statistics 
        % Euclidean distance see formular in document.
        SS1 = ((S_simulated(1) - S_observed(1))/S_simulated(4))^2;
        SS2 = ((S_simulated(2) - S_observed(2))/S_simulated(5))^2;
        SS3 = ((S_simulated(3) - S_observed(3))/S_simulated(6))^2;
        
        % Row 1 of the out vector contains the euclidean distance 
        % the rest of the rows contains the corresponding parameters 
        % used for generating that specific distance.    
        out(1,i) =  ((SS1 + SS2 + SS3)^(0.5));         
        out(2,i) =  param_T;
        out(3,i) =  param_G0;
        out(4,i) =  param_lambda;
        out(5,i) =  param_sigma_N;
        
        disp(i);
        
    end 
    
    % Transpose the out vector in order to use the function "sortrow"
    out = out';
    % Sort the "out" matrix so that the lowest euclidean distance is at the
    % (1,1) matrix position and highest distance is at (max,1) 
    out = sortrows(out); 
    % Transpose matrix back 
    out = out';
    
    % Preallocating vectors
    params_T_lastExtract       = zeros(1,nbr_extract);
    params_G0_lastExtract      = zeros(1,nbr_extract);
    params_lambda_lastExtract  = zeros(1,nbr_extract);
    params_sigma_N_lastExtract = zeros(1,nbr_extract);
    
    params_T       = zeros(1,nbr_extract*iterations);
    params_G0      = zeros(1,nbr_extract*iterations);
    params_lambda  = zeros(1,nbr_extract*iterations);
    params_sigma_N = zeros(1,nbr_extract*iterations);
    
    % Extract the parameters that was used for generating the summary
    % statistics that was closest to the observed summary statistics
    % i.e had the lowest euclidean distances
   
    for i = 1:nbr_extract
        % Following vectors holds the CURRENT extracted parameter values 
        % that was within euclidean distance       
        params_T_lastExtract(i)       = out(2,i);
        params_G0_lastExtract(i)      = out(3,i);
        params_lambda_lastExtract(i)  = out(4,i);
        params_sigma_N_lastExtract(i) = out(5,i);
        
        % Following vectors holds ALL extracted parameter values 
        % that was within euclidean distance
        params_T(index)       = out(2,i);
        params_G0(index)      = out(3,i);
        params_lambda(index)  = out(4,i);
        params_sigma_N(index) = out(5,i);   
        
        index = index + 1; 
    end
    
   % Update the min/max values for the parameters for the next iteration 
   T_a      = min(params_T_lastExtract);
   T_b      = max(params_T_lastExtract);
   G0_a     = min(params_G0_lastExtract);
   G0_b     = max(params_G0_lastExtract);
   lambda_a = min(params_lambda_lastExtract);
   lambda_b = max(params_lambda_lastExtract);
   sigmaN_a = min(params_sigma_N_lastExtract);
   sigmaN_b = max(params_sigma_N_lastExtract);
   
   % Calculate the mean and variance for each parameter from 
   % each iteration, and save in the matrix meanVar_params. 
   % Column 1 contains all the means and variances 
   % from iteration 1, column 2 for iteration 2...
   meanVar_params(1,a) = mean(params_T_lastExtract);
   meanVar_params(2,a) = mean(params_G0_lastExtract);
   meanVar_params(3,a) = mean(params_lambda_lastExtract);
   meanVar_params(4,a) = mean(params_sigma_N_lastExtract);
   meanVar_params(5,a) = var(params_T_lastExtract);
   meanVar_params(6,a) = var(params_G0_lastExtract);
   meanVar_params(7,a) = var(params_lambda_lastExtract);
   meanVar_params(8,a) = var(params_sigma_N_lastExtract);
   % The mean and variances can now be used to plot a normal distribution of the parameters
   % for all iterations.
   
   disp(a);
end 

toc

   
%{    
t = tiledlayout(2,2, 'TileSpacing', 'none', 'Padding', 'compact');
title(t,'Parameter estimation using ABC rejection algorithm for Turin model');

nexttile
scatter(x,params_T)
yline(7.8e-9); % Parameter used for generating observed data 
title('T - reverbation time ');
xlabel('Accepted sample nbr.');

nexttile 
scatter(x,params_G0)
yline(4.07e-9); % Parameter used for generating observed data 
title('G_0 - reverbation gain');
xlabel('Accepted sample nbr.');

nexttile
scatter(x,params_lambda)
yline(10e-12); % Parameter used for generating observed data 
title('\lambda - arrival rate');
xlabel('Accepted sample nbr.');

nexttile
scatter(x,params_sigma_N)
yline(1.673e-4); % Parameter used for generating observed data 
title('\sigma_n - sigma noise');
xlabel('Accepted sample nbr.');

%}


 
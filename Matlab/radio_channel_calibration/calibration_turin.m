clear

%% True data ("measured")
N = 50;% % Number of data sets to generate and average over
B = 4e9; % Bandwidth of signal
Ns = 801; % Number of sample points in each data set
% Parameters
T = 7.8e-9; % Reverberation time
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second
sigma_N = sqrt(0.28e-9); % Nois'e standard deviation

theta_real = [T G0 lambda sigma_N];
S_measured = create_statistics(1,N, theta_real(1), theta_real(2), theta_real(3), theta_real(4), B, Ns);

%% Synthetic likelihood
% Number of summary statistic realisations per iteration
Nr = 20; 
N = 50; % Number of Turin data sets to generate 

% Known parameters
B = 4e9; % Bandwidth of signal
Ns = 801; % Number of sample points in each data set

% Guess of model parameter distribution (uniform) 
Tmax = 10e-9;
Tmin = 5e-9;
G0max = db2pow(-50);
G0min = db2pow(-100);
lambdamax = 12e9;
lambdamin = 5e8;
sigma_Nmax = sqrt(0.25e-10);
sigma_Nmin = sqrt(0.25e-8);

% Make first guess from uniform distribution
T = Tmin + (Tmax - Tmin)*rand; % Reverberation time
G0 = G0min + (G0max - G0min)*rand; % Reverberation gain converted from dB to power
lambda = lambdamin + (lambdamax - lambdamin)*rand; % randomly chosen arrival rate per second lambda 
sigma_N = sigma_Nmin + (sigma_Nmax - sigma_Nmin)*rand; % Noise standard deviation

theta_guess = [T G0 lambda sigma_N]; % Initial guess

% Generating statistics of synthetic data using initial guess on theta
S_star = create_statistics(Nr, N, theta_guess(1), theta_guess(2), theta_guess(3), theta_guess(4), B, Ns);

% S_star_scaled = (S_star - mean(S_star,1))./std(S_star,1); % Z-score normalisation
% S_measured_scaled = (S_measured - mean(S_star,1))./std(S_star,1);
%% Tester
% mu = mean(S_star); % Calculate mean of summary statistics from simulated data 
% Sigma = cov(S_star); % Calculate covariance of summary statistics from simulated data 
% L = -1/2*( ( (S_measured-mu)/Sigma) ) * (S_measured-mu)' - 1/2*log(det(Sigma)); % Synthetic likelihood L

%% Generate first likelihood based on initial guess

L_guess = synth_loglikelihood(S_measured, S_star);
M = 100; % Number of BSL posterior samples 
theta = zeros(M,4); % Buffer for parameter values
L = zeros(M,1); % Buffer for log likelihood values
theta(1,:) = theta_guess; % First guess
L(1) = L_guess; % First log likelihood

%% Monte Carlo
for i = 2:M
    % Propose new parameter set theta + q, q = uniform
    T = Tmin + (Tmax - Tmin)*rand; % Reverberation time
    G0 = G0min + (G0max - G0min)*rand; % Reverberation gain converted from dB to power
    lambda = lambdamin + (lambdamax - lambdamin)*rand; % randomly chosen arrival rate per second lambda 
    sigma_N = sigma_Nmin + (sigma_Nmax - sigma_Nmin)*rand; % Noise standard deviation
    theta_prop = [T G0 lambda sigma_N];
    
    % Statistics of Nr data sets with proposed parameter set
    S_star = create_statistics(Nr, N, theta_prop(1), theta_prop(2), theta_prop(3), theta_prop(4), B, Ns);
    % Calculate new log likelihood 
    L_prop = synth_loglikelihood(S_measured, S_star);
    % Compare with old log likelihood
    if (exp(L_prop - L_guess) > rand) % Metropolis Hastings acceptance criterion (P_theta*/P_theta)
        theta_guess = theta_prop;
        L_guess = L_prop;
    end
    % Save parameter 
    theta(i,:) = theta_guess;
    L(i) = L_guess;
end


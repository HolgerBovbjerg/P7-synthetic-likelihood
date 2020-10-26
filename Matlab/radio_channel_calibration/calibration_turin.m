clear


%% True data ("measured")
N = 1;% % Number of data sets to generate and average over
B = 4e9; % Bandwidth of signal: 4 GHz
Ns = 801; % Number of sample points in each data set
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second
sigma_N = sqrt(0.28e-9); % Noise standard deviation

theta_real = [T G0 lambda sigma_N];

S_measured = mean(create_statistics(N, T, G0, lambda, sigma_N, B, Ns),1);


%% Synthetic likelihood
% Known parameters
N = 200;% % Number of data sets to generate and average over
B = 4e9; % Bandwidth of signal
Ns = 801; % Number of sample points in each data set

% Guess of model parameters theta 
Tmax = 9e-9;
Tmin = 7e-9;
T = Tmax + (Tmax - Tmin)*rand; % Reverberation time
G0max = db2pow(-70);
G0min = db2pow(-90);
G0 = G0max + (G0max - G0min)*rand; % Reverberation gain converted from dB to power
lambdamax = 12e9;
lambdamin = 8e9;
lambda = lambdamax + (lambdamax - lambdamin)*rand; % randomly chosen arrival rate per second lambda 
sigma_Nmax = sqrt(0.35e-9);
sigma_Nmin = sqrt(0.25e-9);
sigma_N = sigma_Nmax + (sigma_Nmax - sigma_Nmin)*rand; % Noise standard deviation

theta_guess = [T G0 lambda sigma_N];

% Generating synthetic data using initial guess on theta
S = create_statistics(N, theta_guess(1), theta_guess(2), theta_guess(3), theta_guess(4), B, Ns);

S_star = mean(S,1);

L = synth_loglikelihood(S_measured, S_star);

% from L update theta using MCMC?

% calcuate untill L is maximised given some epsilon or some number of
% iterations



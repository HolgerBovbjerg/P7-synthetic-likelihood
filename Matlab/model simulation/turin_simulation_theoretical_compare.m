%% clear
clear 

%% Initial choices made on model. 
% We have chosen some values based on "Estimator for Stochastic Channel Model without
% Multipath Extraction using Temporal Moments" by Ayush Bharti et al. 

N = 200;% % Number of data sets to generate and average over
B = 4e9; % Bandwidth of signal: 4 GHz
Ns = 801; % Number of sample points in each data set
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second
sigma_N = sqrt(0.28e-9); % Noise standard deviation

[P_y, Y_k, t] = sim_turin_matrix_gpu(N,B,Ns,T,G0,lambda,sigma_N);

% We use the formulas from the paper before. 
P_h_theoretical = G0*exp(-(t/T));

% We use P_Y = E_s * P_h + noise (Noise is already included in simulation)
P_y_simulated = P_y; %+ sigma_N^2/Ns;
P_y_theoretical = P_h_theoretical + sigma_N^2/Ns; % theoretical does not need bandwidth scaling


%% Generation of plots showing the power spectrum. 
figure
plot(t*1e9,pow2db(P_y_theoretical), 'DisplayName', "P_y theoretical")
hold on
plot(t*1e9,pow2db(P_y_simulated), 'DisplayName', "P_y simulated")
xlim([0 200])
ylim([-130 -80])
xlabel("Time [ns]")
ylabel("Power [dB]")
lgd = legend;

exportgraphics(gcf,'turin_theoretical_simulation_compare.pdf','ContentType','vector')
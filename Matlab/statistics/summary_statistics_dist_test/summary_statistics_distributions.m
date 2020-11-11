%% Distribution of summary statistics
clear
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
sigma_N = sqrt(28e-9); % Noise variance
B = 4e9; % Bandwidth of signal: 4 GHz
Ns = 801; % Number of frequency samples in transfer function

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second

Nr = 1000;
S = zeros(4,Nr);
tic
temp = create_statistics(10, T, G0, lambda, sigma_N, B, Ns);
t_per_run = toc/10;


est_time = Nr*t_per_run

%%
tic
S = create_statistics(Nr, T, G0, lambda, sigma_N, B, Ns);
toc
%%
bins = 15;
figure(1)
plotstatistics =  tiledlayout(1,4);
title(plotstatistics,"Histograms of summary statistics, from "+ Nr +" summary statistic vectors") 

nexttile
hist(S(:,1),bins)
title("Histogram of mean of 0th temporal moment")

nexttile
hist(S(:,2),bins)
title("Histogram of mean of 1st temporal moment")

nexttile
hist(S(:,3),bins)
title("Histogram of mean of 2nd temporal moment")

nexttile
hist(S(:,4),bins)
title("Histogram of mean of 3rd temporal moment")

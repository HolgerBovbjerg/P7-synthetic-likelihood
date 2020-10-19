%% Distribution of summary statistics
clear
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
sigma_N = sqrt(28e-9); % Noise variance
B = 4e9; % Bandwidth of signal: 4 GHz
Ns = 801; % Number of frequency samples in transfer function

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second

Nr = 50;    % Number of Turin realizations pr summary statistic vector
Nl = 1000;   % Number of different summary statistic vectors generated

ss = zeros(Nl,8); % Matrix for summary statistics

tic
temp = create_statistics(Nr, T, G0, lambda, sigma_N, B, Ns);
t_per_stat = toc/Nr;


est_time = Nl*t_per_stat

%%

T_min = 7.8e-11; T_max = 7.8e-7;
T = T_min + (T_max-T_min).*rand(Nl,1);

G0_min = db2pow(-110); G0_max = db2pow(-60);
G0 = G0_min + (G0_max-G0_min).*rand(Nl,1);

sigma_N_min = sqrt(5e-10); sigma_N_max = sqrt(5e-9);
sigma_N = sigma_N_min + (sigma_N_max-sigma_N_min).*rand(Nl,1);

lambda_min = 50e6; lambda_max = 50e7;
lambda = lambda_min + (lambda_max-lambda_min).*rand(Nl,1);

%%
tic
for i = 1:Nl    
    ss(i,:) = create_statistics(Nr, T(i), G0(i), lambda(i), sigma_N(i), B, Ns);
end
toc
%%
bins = 85;
figure(1)
plotstatistics =  tiledlayout(1,5);
title(plotstatistics,"Histograms of summary statistics, from "+ Nl +" summary statistic vectors") 

nexttile
hist((ss(:,1)),bins)
title("Histogram of mean of 0th temporal moment")

nexttile
hist(ss(:,2),bins)
title("Histogram of mean of 1st temporal moment")

nexttile
hist(ss(:,3),bins)
title("Histogram of mean of 2nd temporal moment")

nexttile
hist(ss(:,4),bins)
title("Histogram of mean of 3rd temporal moment")

nexttile
hist(ss(:,5),bins)
title("Histogram of variance of 0th temporal moment")

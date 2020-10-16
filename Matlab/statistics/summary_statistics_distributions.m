%% Distribution of summary statistics
clear
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
sigma_N = sqrt(28e-9); % Noise variance
N = 10; 
B = 4e9; % Bandwidth of signal: 4 GHz
Ns = 801;
% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second
Nr = 10;
S = zeros(8,Nr);
t_per_run = 0.005;
est_time = N*Nr*t_per_run;

%%
temporal0 = gpuArray(zeros(N,Nr));
temporal1 = gpuArray(zeros(N,Nr));
temporal2 = gpuArray(zeros(N,Nr));
temporal3 = gpuArray(zeros(N,Nr));

for k = 1:Nr
    
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, T, G0, lambda, sigma_N);
    
    t = gpuArray(t);
    temporal0(:,k) = trapz(t,(t.^0.*Pv));
    temporal1(:,k) = trapz(t,(t.^1.*Pv));
    temporal2(:,k) = trapz(t,(t.^2.*Pv));
    temporal3(:,k) = trapz(t,(t.^3.*Pv));
    
    S(:,k) = [mean(temporal0(:,k)); ...
        mean(temporal1(:,k)); ...
        mean(temporal2(:,k)); ...
        mean(temporal3(:,k)); ...
        var(temporal0(:,k)); ...
        var(temporal1(:,k)); ...
        var(temporal2(:,k)); ...
        var(temporal3(:,k)) ];      
end
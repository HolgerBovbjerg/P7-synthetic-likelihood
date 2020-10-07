%% clear
clear 

%% 

%% Choices made on model
N = 625; % Number of data sets to generate
B = 4e9; % Bandwidth of signal
Ns = 801; % Number of sample points
deltaf = B/(Ns-1); % Frequency seperation
tmax = 1/deltaf; 
T = 7.8e-9; % Reverberation time
G0 = db2pow(-83.9); % Reverberation gain

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10; % randomly chosen arrival rate lambda
ldist = makedist('poisson',lambda);


%% Simulate model
Hk = zeros(Ns,N); % buffer for generated channel response data
sigma_N = sqrt(28e-9); % Noise variance

for n = 1:N
    % Ns = random(kdist,1,1)
    lmax = random(ldist,1,1); % Number of multipath components
    tau = rand(lmax,1)*tmax; % time-delays
    tau = sort(tau);
    for l = 1:length(lmax) % For every multipath component
        sigma_alpha = sqrt(G0*exp(-(tau(l)/T) ) / lambda); % Calculate variance
        alphadist = makedist('rayleigh', sigma_alpha);
        phidist = makedist('uniform', 0, 2*pi);
        alpha = random(alphadist,1,1);
        phi = random(phidist,1,1);
        for k = 1:Ns % For every frequency index
            Hk(k,n) = Hk(k,n) + alpha*exp(-1j*phi)*exp(-1j*2*pi*deltaf*k*tau(l));
        end
    end
end

%% Mean over N realizations
Hkmean = zeros(Ns,1);
for i = 1:Ns
    Hkmean(i) = mean(Hk(i,:));
end

%% Calculate Power 
P_h_simulated = abs(ifft(Hkmean)).^2;

t = (0:Ns-1)./(deltaf*Ns);

P_h_theoretical = G0*exp(-(t/T));

% P_y_simulated = B*P_h_simulated + sigma_N^2/Ns;

%% Plots
tiledlayout(2,1)
nexttile
plot(t*1e9,pow2db(P_h_theoretical))
xlim([0 200])
xlabel("Time [ns]")
% ylabel("Power [dB]")
nexttile
plot(t*1e9,pow2db(P_h_simulated))
xlim([0 200])
xlabel("Time [ns]")
ylabel("Power [dB]")

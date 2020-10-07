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
lambda = 100; % randomly chosen arrival rate lambda
ldist = makedist('poisson',lambda);


%% Simulate model
Hk = zeros(Ns,N); % buffer for generated channel response data
sigma_N = sqrt(28e-9); % Noise variance

for n = 1:N
    % Ns = random(kdist,1,1)
    lmax = random(ldist,1,1); % Number of multipath components
    tau = rand(lmax,1)*tmax; % time-delays
    tau = sort(tau);
    for k = 1:Ns % For every frequency index
        for l = 1:length(lmax) % For every multipath component
            sigma_alpha = sqrt(G0*exp(-(tau(l)/T) ) / lambda); % Calculate variance
            alphadist = makedist('rayleigh', sigma_alpha);
            alpha_real = random(alphadist,1,1);
            alpha_imag = random(alphadist,1,1);
            alpha = 1/sqrt(2)*(alpha_real+1j*alpha_imag);
            Hk(k,n) = Hk(k,n) + alpha*exp(-1j*2*pi*deltaf*k*tau(l));
        end
    end
end

%%
Hkmean = zeros(Ns,1);
for i = 1:Ns
    Hkmean(i) = mean(Hk(i,:));
end

% P_h_theoretical = mean(abs(ifft(Hk)).^2,2);
%%
P_h_simulated = abs(ifft(Hkmean)).^2;

t = (0:Ns-1)./(deltaf*Ns);

P_h_theoretical = G0*exp(-(t/T));

tiledlayout(2,1)
nexttile
plot(t*1e9,P_h_theoretical)
xlim([0 200])
xlabel("Time [ns]")
% ylabel("Power [dB]")
nexttile
plot(t*1e9,pow2db(P_h_simulated))
xlim([0 200])
xlabel("Time [ns]")
ylabel("Power [dB]")
%% Compute ensemble mean of amplitude at delay t over N realisations
% rhomean = zeros(100,1); % buffer for ensemble mean of channel response
% for t=1:length(rhomean) %for every time delay 
%     for n = 1:N % for every experiment
%         rhomean(t) = rhomean(t) + rho(t,n); % sum channel response
%     end
%     rhomean(t) = rhomean(t)/N; % Divide by number of experiments
% end
%         
%% Plotting data

% trange = (0:99); % range of time delays  
% 
% plot(trange,abs(rhomean))
% title("Simulated response of radio channel using Turin model")
% xlabel("Time delay (arbitrary)")
% ylabel("Mean amplitude over " +  N + " samples with " + kmax + " rays")
% % xlim([lambda-lambda*3/4 lambda+lambda*3/4])

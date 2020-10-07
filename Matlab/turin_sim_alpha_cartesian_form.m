%% clear
clear 

%% 

%% Initial choices made on model. 
% We have chosen some values based on "Estimator for Stochastic Channel Model without
% Multipath Extraction using Temporal Moments" by Ayush Bharti et al. 

N = 625; % Number of data sets to generate and average over
B = 4e9; % Bandwidth of signal: 4 GHz
Ns = 801; % Number of sample points in each data set
deltaf = B/(Ns-1); % Frequency seperation: 5 MHz
tmax = 1/deltaf; % Maximum delay, found from the bandwidth via frequency seperation
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 100; % randomly chosen arrival rate lambda
ldist = makedist('poisson',lambda); 


%% Simulate model
% We then simulate the transfer function, H_k, of the channel, using the
% Turin model. 

Hk = zeros(Ns,N); % buffer for generated channel response data
sigma_N = sqrt(28e-9); % Noise variance

% We run the simulation N times, creating new data sets for each
% realization. 
for n = 1:N
    % Ns = random(kdist,1,1) 
    lmax = random(ldist,1,1);   % Number of multipath components created from the Poisson distribution.
    tau = rand(lmax,1)*tmax;    % time-delays drawn uniformly between 0 and the maximum delay.  
    tau = sort(tau);            % The time delays are sorted, for easier plotting later. 
    alpha = zeros(lmax,1);      % Buffer for generated alpha values.
   
    % For every multipath component a complex gain is generated, based on a
    % sigma generated from a corresponding delay time value. 
    % The complex number is generated in cartesian form by drawing the real
    % and the imaginary part seperately from a Rayleigh distribution. 
    for i = 1:length(lmax) 
            sigma_alpha = sqrt(G0*exp(-(tau(i)/T) ) / lambda); % Calculate variance
%             alphadist = makedist('normal',0, sigma_alpha);
            alphadist = makedist('rayleigh', 0, sigma_alpha);
            alpha_real = random(alphadist,1,1);
            alpha_imag = random(alphadist,1,1);
            alpha(i) = 1/sqrt(2)*(alpha_real+1j*alpha_imag);           
    end
    
    % For every frequency index, k, the transfer function is calculated as
    % a sum with contributions from each multipath component. 
    for k = 1:Ns 
        for l = 1:length(lmax) % For every multipath component
            Hk(k,n) = Hk(k,n) + alpha(l)*exp(-1j*2*pi*deltaf*k*tau(l));
        end
    end
end

%% Averaging over the N realizations

Hkmean = zeros(Ns,1);
for i = 1:Ns
    Hkmean(i) = mean(Hk(i,:));
end


%% Simulating the power spectrum, P_h, from the transfer function
% We use the formulas from the paper before. 

P_h_simulated = abs(ifft(Hkmean)).^2;
P_h_theoretical = G0*exp(-(t/T));


% Generate times for time axis in plot
t = (0:Ns-1)./(deltaf*Ns);



% Generation of plots showing the power spectrum. 
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

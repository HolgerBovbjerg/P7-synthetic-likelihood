function P_y_simulated = turin_sim_alpha_cartesian_form_claus(lambda,G0,T,N)

%% clear
%%clear 

%% Initial choices made on model. 
% We have chosen some values based on "Estimator for Stochastic Channel Model without
% Multipath Extraction using Temporal Moments" by Ayush Bharti et al. 

%N = 625; % Number of data sets to generate and average over
B = 4e9; % Bandwidth of signal: 4 GHz
Ns = 801; % Number of sample points in each data set
deltaf = B/(Ns-1); % Frequency seperation: 5 MHz
tmax = 1/deltaf; % Maximum delay, found from the bandwidth via frequency seperation
%T = 7.8e-9; % Reverberation time: 7.8 ns
%G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
sigma_N = sqrt(28e-9); % Noise variance

% Time delay tau is a possion arrival process with mean delay lambda
%lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second

%% Simulate model
% We then simulate the transfer function, H_k, of the channel, using the
% Turin model. 

Hk = zeros(Ns,N); % buffer for generated channel response data
P_h = zeros(Ns,N); % Buffer for Power of impulse response
ldist = makedist('poisson',tmax*lambda); % distribution for number 
                                             % of multipaths to include is 
                                             % a function of lambda and
                                             % maximum delay
% We run the simulation N times, creating new data sets for each
% realization. 
for n = 1:N
    lmax = random(ldist,1,1);   % Number of multipath components, created from the Poisson distribution.
    tau = rand(lmax,1)*tmax;    % time-delays, drawn uniformly between 0 and the maximum delay.  
    tau = sort(tau);
    % For every multipath component a complex gain is generated, based on a
    % sigma generated from a corresponding delay time value. 
    % The complex number is generated in cartesian form by drawing the real
    % and the imaginary part seperately from a normal distribution. 
    for i = 1:lmax 
            sigma_alpha = sqrt(G0*exp(-(tau(i)/T) ) / lambda);% Calculate variance using eq 13 and Ph = Lambda*sigma^2 from Ayush paper.
            
            %% 1
            alphadist = makedist('normal',0,sigma_alpha);
            alpha_real = random(alphadist,1,1);
            alpha_imag = random(alphadist,1,1);
            % The complex valued alpha is created by combining the real and
            % imaginary parts. 
            alpha = 1/sqrt(2)*(alpha_real+1j*alpha_imag);  
            % For every frequency index, k, the contribution from multipath
            % component, i, is added to the transfer function. 
            for k = 1:Ns 
                Hk(k,n) = Hk(k,n) + alpha*exp(-1j*2*pi*deltaf*k*tau(i));
            end
            
    end
    
    P_h(:,n) = abs(ifft(Hk(:,n))).^2;
end
%% Averaging over the N realizations

P_h_mean = mean(P_h,2);


% P_h_mean = zeros(Ns,1);
% for i = 1:Ns
%     P_h_mean(i) = mean(P_h(i,:));
% end

%% Simulating the power spectrum, P_h, from the transfer function
% Generate timestamps, in seconds, for time axis in plot
t = (0:Ns-1)./(deltaf*Ns);

% We use the formulas from the paper before. 
P_h_simulated = P_h_mean;
P_h_theoretical = G0*exp(-(t/T));

% We use P_Y = E_s * P_h + noise
P_y_simulated = B*P_h_simulated + sigma_N^2/Ns;
P_y_theoretical = P_h_theoretical + sigma_N^2/Ns;

% Generation of plots showing the power spectrum. 
plot(t*1e9,pow2db(P_y_theoretical), 'DisplayName', "P_y theoretical")
hold on
plot(t*1e9,pow2db(P_y_simulated), 'DisplayName', "P_y simulated")
title("P_y simulated")
xlim([0 200])
ylim([-110 -80])
xlabel("Time [ns]")
ylabel("Power [dB]")
lgd = legend;


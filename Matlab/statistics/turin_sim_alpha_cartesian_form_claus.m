function [P_y_simulated P_h alpha tau, t] = turin_sim_alpha_cartesian_form_claus(lambda,G0,T,N,sigma_N)

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
%sigma_N = sqrt(28e-9); % Noise variance

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
lmax = poissrnd(tmax*lambda);   % Number of multipath components, created from the Poisson distribution.
alpha = zeros(N,lmax);
tau = zeros(lmax,N);   
for n = 1:N

        tau(:,n) = rand(lmax,1)*tmax;    % time-delays, drawn uniformly between 0 and the maximum delay.  
        % For every multipath component a complex gain is generated, based on a
        % sigma generated from a corresponding delay time value. 
        % The complex number is generated in cartesian form by drawing the real
        % and the imaginary part seperately from a normal distribution. 
        sigma_alpha = sqrt(G0*exp(-(tau(:,n)*(1/T)) ) / lambda);% Calculate variance using eq 13 and Ph = Lambda*sigma^2 from Ayush paper.
        % The complex valued alpha is created by combining the real and
        % imaginary parts. 
        alpha(n,:) = sigma_alpha * 1/sqrt(2) .* (randn(lmax,1) +  1j*randn(lmax,1));  
        % For every frequency index, k, the contribution from every multipath
        % component is added to the transfer function. 
        k = (1:Ns);
        Hk(:,n) = (exp(-1j*2*pi*deltaf*k.*tau(:,n)).' * alpha(n,:)');
    end
%% Averaging over the N realizations
   
% Generate noise vector
noise = sigma_N^2 * (1/sqrt(2))* (randn(Ns,N) + 1j*randn(Ns,N));
% Add noise to transfer function with complex normal dist
Hk = Hk + noise;
% Power delay profile
P_h = abs(ifft(Hk,[],1)).^2;
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
P_y_simulated = B*P_h_simulated;
P_y_theoretical = P_h_theoretical + sigma_N^2/Ns;


% Generation of plots showing the power spectrum. 
figure(1)
% plot(t*1e9,pow2db(P_y_theoretical), 'DisplayName', "P_y theoretical")
 hold on
 plot(t*1e9,pow2db(P_y_simulated), 'DisplayName', "P_y simulated")
 title("P_y simulated")
 %xlim([0 200])
% ylim([-110 -80])
 xlabel("Time [ns]")
 ylabel("Power [dB]")
 lgd = legend;


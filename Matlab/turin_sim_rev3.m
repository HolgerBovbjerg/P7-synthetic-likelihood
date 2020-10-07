%% clear
clear 
%%
kmax = 30; % Number of rays to include
N = 10; % Number of data sets to generate

% Phase phi is a uniform r.v. between 0 and 2*pi
phidist = makedist('uniform',0,2*pi);

% Gain a is a lognormal r.v. 
mu = 0.6; % randomly chosen mu
sigma = 1; % randomly chosen sigma
adist = makedist('Lognormal', mu, sigma);
% adist = makedist('Rayleigh', sigma);

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 50; % randomly chosen lambda
taudist = makedist('Exponential',1/lambda);


rho = zeros(kmax,N); % buffer for generated channel response data
time = zeros(kmax,N);
for n = 1:N
    phi = random(phidist,kmax,1); % Generate phase vector
    a = random(adist,kmax,1); % Generate amplitude vector
    time(:,n) = random(taudist,kmax,1); % Generate delay vector
    time = sort(time);
    for k = 1:kmax
        rho(k,n) = a(k) * exp(1j*phi(k)); % Add ray contribution to sum
    end
end




%% Compute ensemble mean of amplitude at delay t over N realisations
rhomean = zeros(kmax,1); % buffer for ensemble mean of channel response
for t=1:length(rhomean) %for every time delay 
    for n = 1:N % for every experiment
        rhomean(t) = rhomean(t) + rho(t,n); % sum channel response
    end
    rhomean(t) = rhomean(t)/N; % Divide by number of experiments
end
        
%% Plotting data

plot(time,abs(rho))
title("Simulated response of radio channel using Turin model")
xlabel("Time delay (arbitrary)")
ylabel("Mean amplitude over " +  N + " samples with " + kmax + " rays")
% xlim([lambda-lambda*3/4 lambda+lambda*3/4])

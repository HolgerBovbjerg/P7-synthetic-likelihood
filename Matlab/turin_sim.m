%% clear
clear 
%%
kmax = 100000; % Number of rays to include
N = 100; % Number of data sets to generate

% Phase phi is a uniform r.v. between 0 and 2*pi
phidist = makedist('uniform',0,2*pi);

% Gain a is a lognormal r.v. 
mu = 0; % randomly chosen mu
sigma = 1; % randomly chosen sigma
adist = makedist('Lognormal', mu, sigma);
% adist = makedist('Rayleigh', sigma);

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10; % randomly chosen lambda
taudist = makedist('poisson',lambda);


rho = zeros(100,N); % buffer for generated channel response data
for n = 1:N
    phi = random(phidist,kmax,1); % Generate phase vector
    a = random(adist,kmax,1); % Generate amplitude vector
    tau = random(taudist,kmax,1); % Generate delay vector
    for t = (1+min(tau)):(1+max(tau)) % For every time delay
        for k = 1:kmax % For every ray "k"
            if (t-1-tau(k) == 0) % if the time delay of k'th ray == current time delay
                rho(t,n) = rho(t,n) + a(k) * exp(1j*phi(k)); 
            end
        end
    end
end


%% Compute ensemble mean of amplitude at delay t over N realisations
rhomean = zeros(100,1); % buffer for ensemble mean of channel response
for t=1:length(rhomean) %for every time delay 
    for n = 1:N % for every experiment
        rhomean(t) = rhomean(t) + rho(t,n); % sum channel response
    end
    rhomean(t) = rhomean(t)/N; % Divide by number of experiments
end
        
%% Plotting data

trange = (0:99); % range of time delays  

plot(trange,abs(rhomean))
title("Simulated response of radio channel using Turin model")
xlabel("Time delay (arbitrary)")
ylabel("Mean amplitude over " +  N + " samples with " + kmax + " rays")
% xlim([lambda-lambda*3/4 lambda+lambda*3/4])

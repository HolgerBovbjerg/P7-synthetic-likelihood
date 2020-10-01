%% clear
clear 
%%
kmax = 1000; % Number of rays to include
N = 100; % Number of data sets to generate

% Generating phase phi as uniform r.v. between 0 and 2*pi
phidist = makedist('uniform',0,2*pi);


% Generating gain a as lognormal r.v. 
mu = 0; % randomly chosen mu
sigma = 1; % randomly chosen sigma
adist = makedist('Lognormal', mu, sigma);
% adist = makedist('Rayleigh', sigma);


% Generate t as possion arrival process with mean delay lambda
lambda = 10; % randomly chosen lambda
tdist = makedist('poisson',lambda);


rho = zeros(kmax,N);
for n = 1:N
    phi = random(phidist,kmax,1); %phase of k'th ray
    a = random(adist,kmax,1);
    tau = random(tdist,kmax,1);
    for t = (1+min(tau)):max(tau)
        for k = 1:kmax
            if (t-1-tau(k) == 0)
                rho(t,n) = rho(t,n) + a(k) + exp(1j*phi(k));
            end
        end
    end
end


%% Compute ensemble mean of amplitude at delay t over N realisations
rhomean = zeros(kmax,1);
for t=1:kmax
    for n = 1:N
        rhomean(t) = rhomean(t) + rho(t,n);
    end
    rhomean(t) = rhomean(t)/N;
end
        
trange = (0:kmax-1); 
plot(trange,abs(rhomean))
title("Simulated response of radio channel using Turin model")
xlabel("Time delay (arbitrary)")
ylabel("Mean amplitude over " +  N + " samples with " + kmax + " rays")
xlim([lambda-lambda*3/4 lambda+lambda*3/4])

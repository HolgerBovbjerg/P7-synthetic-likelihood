%% clear
clear 
%%
tMax = 15; % Maximum time delay
N = 100; % Number of data sets to generate

% Phase phi is a uniform r.v. between 0 and 2*pi
phidist = makedist('uniform',0,2*pi);

% Gain a is a lognormal r.v. 
mu = 0.6; % randomly chosen mu
sigma = 1; % randomly chosen sigma
adist = makedist('Lognormal', mu, sigma);
% adist = makedist('Rayleigh', sigma);

% Time delay tau is a possion arrival process with mean arival rate lambda
lambda = 50; % randomly chosen lambda
taudist = makedist('poisson',lambda);
% taudist = makedist('Exponential',lambda);

time = 0;

for t = 1:tMax % For every time delay
    numPoints = random(taudist,1,1);
    time = [time; rand(numPoints,1) + (t-1)];
end

time = sort(time);

tdelay=zeros(length(time),1);
tdelay(1) = time(1);
for t = 2:length(time)
    tdelay(t) = time(t) - time(t-1);
end

a = random(adist,length(time),1);
phi = random(phidist,length(time),1);

rho = a.*exp(1j*phi);

fitdist(tdelay,'exponential')
plot(time,abs(rho))



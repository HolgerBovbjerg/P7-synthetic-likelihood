%% clear
clear

%% Generate tau
% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10; % randomly chosen lambda
taudist = makedist('Exponential',lambda);
tau = random(taudist,10000,1);

histogram(tau)
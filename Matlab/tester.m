%% clear
clear 
%%
kmax = 15; % Number of rays to include
N = 1; % Number of data sets to generate

% Phase phi is a uniform r.v. between 0 and 2*pi
phidist = makedist('uniform',0,2*pi);

% Gain a is a lognormal r.v. 
mu = 0; % randomly chosen mu
sigma = 1; % randomly chosen sigma
adist = makedist('Lognormal', mu, sigma);
% adist = makedist('Rayleigh', sigma);

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 4; % randomly chosen lambda
% taudist = makedist('poisson',lambda);
taudist = makedist('Exponential',lambda);

taumax = 300;
timestep = 1;

rho = zeros(taumax/timestep,N); % buffer for generated channel response data

for n = 1:N
    phi = random(phidist,kmax,1); % Generate phase vector
    a = random(adist,kmax,1); % Generate amplitude vector
    tau = random(taudist,kmax,1); % Generate delay vector
    for t = 0:timestep:taumax % For every time delay
        for k = 1:kmax % For every ray "k"
            if (tau(k) <= t ) && ( (t-timestep) < t )  % if the time delay of k'th ray == current time delay
                rho(t/timestep) = rho(t/timestep) + a(k) * exp(1j*phi(k)); % Add ray contribution to sum
            end
        end
    end
end

rhomean = zeros(taumax,1); % buffer for ensemble mean of channel response
t = (0:timestep:taumax-1);

plot(t, abs(rho))

% 
% for t=1:length(rhomean) %for every time delay 
%     for n = 1:N % for every experiment
%         rhomean(t) = rhomean(t) + rho(t,n); % sum channel response
%     end
%     rhomean(t) = rhomean(t)/N; % Divide by number of experiments
% end
%        
% plot(abs(rhomean))
% xlim([0 60])
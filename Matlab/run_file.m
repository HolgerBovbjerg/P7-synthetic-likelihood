T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second

n = 50;


turin_sim_alpha_cartesian_form(lambda,G0,T,n)

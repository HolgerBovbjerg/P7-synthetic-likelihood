clear all
tic
%T = 7.8e-9;
T = linspace(5e-9, 10e-9, 50); % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second

n = 5;
mean0 = zeros(1,50);
for i = 1:length(T)
    P = (turin_sim_alpha_cartesian_form_claus(lambda,G0,T(i),n));
    mean0(i) = mean(P);
    i
end
toc
figure(2)
plot(T,mean0,"o")

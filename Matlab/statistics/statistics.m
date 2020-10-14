% Matlab script for finding summary statistics from samples from the Turin
% model. In eact section one parameter is varied while the others are held
% fixed. 50 different values of each parameter is simulated

% The temporal moments are found using equation 3 in the paper "Parameter
% estimation for stochastic channel models using temporal moments"

% When varying a value, n iterations are being done. The mean and variance
% of these iterations are then found as our summary statistics.

%% Varying T
clear all

tic % Starting a timer to see how long it takes to simulate
T = linspace(7.8e-11,7.8e-9,50); % Reverberation time, 
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
n = 300; % Number of iterations 

sigma_noise = sqrt(0.28e-9); % noise variance
order = 3;
for k = 1:length(T)
    [P Pv alpha tau] = turin_sim_alpha_cartesian_form_claus(lambda,G0,T(k),n,sigma_noise);
    for l = 1:n
            %   0th order moment
        moment0bef(l) = abs(alpha(l,:)).^2*tau(:,l).^0; % Vector containing the moment of each iteration
        moment0(k) = mean(moment0bef);                  % Finding the mean of the moment of the n iterations
        var0(k) = var(moment0bef);                      % Finding the variance of the moment of the n iterations
            %   1st order moment
        moment1bef(l) = abs(alpha(l,:)).^2*tau(:,l).^1;
        moment1(k) = mean(moment1bef);
        var1(k) = var(moment1bef);
            %   2nd order moment
        moment2bef(l) = abs(alpha(l,:)).^2*tau(:,l).^2;
        moment2(k) = mean(moment2bef);
        var2(k) = var(moment2bef);
            %   3rd order moment
        moment3bef(l) = abs(alpha(l,:)).^2*tau(:,l).^3;
        moment3(k) = mean(moment3bef);
        var3(k) = var(moment3bef);
    end
    k
end    
toc
figure(2)
plotT = tiledlayout(4,2)
title(plotT,"Summary statistics varying T")
nexttile
plot(T,moment0,"o")
title("Mean (m0)")

nexttile
plot(T,var0,"o")
title("Var (m0)")

nexttile
plot(T,moment1,"o");
title("Mean (m1)")

nexttile
plot(T,var1,"o")
title("Var (m1)")

nexttile
plot(T,moment2,"o");
title("Mean (m2)")

nexttile
plot(T,var2,"o")
title("Var (m2)")

nexttile
plot(T,moment3,"o");
title("Mean (m3)")

nexttile
plot(T,var3,"o")
title("Var (m2)")



%% Varying G0
clear all

tic
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = linspace(db2pow(-110), db2pow(-60), 50); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
n = 300;
sigma_noise = sqrt(0.28e-9);
order = 3;
for k = 1:length(G0)
    [P Pv alpha tau] = turin_sim_alpha_cartesian_form_claus(lambda,G0(k),T,n,sigma_noise);
    for l = 1:n
        moment0bef(l) = abs(alpha(l,:)).^2*tau(:,l).^0;
        moment0(k) = mean(moment0bef);
        var0(k) = var(moment0bef);
        
        moment1bef(l) = abs(alpha(l,:)).^2*tau(:,l).^1;
        moment1(k) = mean(moment1bef);
        var1(k) = var(moment1bef);
        
        moment2bef(l) = abs(alpha(l,:)).^2*tau(:,l).^2;
        moment2(k) = mean(moment2bef);
        var2(k) = var(moment2bef);
        
        moment3bef(l) = abs(alpha(l,:)).^2*tau(:,l).^3;
        moment3(k) = mean(moment3bef);
        var3(k) = var(moment3bef);
    end
    k
end    
toc
figure(3)
plotT = tiledlayout(4,2)
title(plotT,"Summary statistics varying G0")

nexttile
plot(G0,moment0,"o")
title("Mean (m0)")

nexttile
plot(G0,var0,"o")
title("Var (m0)")

nexttile
plot(G0,moment1,"o");
title("Mean (m1)")

nexttile
plot(G0,var1,"o")
title("Var (m1)")

nexttile
plot(G0,moment2,"o");
title("Mean (m2)")

nexttile
plot(G0,var2,"o")
title("Var (m2)")

nexttile
plot(G0,moment3,"o");
title("Mean (m3)")

nexttile
plot(G0,var3,"o")
title("Var (m3)")


%% Varying noise variance
clear all

tic
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
n = 30;
sigma_noise = linspace(sqrt(0.28e-11),sqrt(0.28e-8),50);
order = 3;

for k = 1:length(sigma_noise)
    [P Pv alpha tau] = turin_sim_alpha_cartesian_form_claus(lambda,G0,T,n,sigma_noise(k));
    for l = 1:n
        moment0bef(l) = abs(alpha(l,:)).^2*tau(:,l).^0;
        moment0(k) = mean(moment0bef);
        var0(k) = var(moment0bef);
        
        moment1bef(l) = abs(alpha(l,:)).^2*tau(:,l).^1;
        moment1(k) = mean(moment1bef);
        var1(k) = var(moment1bef);
        
        moment2bef(l) = abs(alpha(l,:)).^2*tau(:,l).^2;
        moment2(k) = mean(moment2bef);
        var2(k) = var(moment2bef);
        
        moment3bef(l) = abs(alpha(l,:)).^2*tau(:,l).^3;
        moment3(k) = mean(moment3bef);
        var3(k) = var(moment3bef);
    end
    k
end    
toc
figure(4)
plotT = tiledlayout(4,2)
title(plotT,"Summary statistics varying noise variance")

nexttile
plot(sigma_noise,moment0,"o")
title("Mean (m0)")

nexttile
plot(sigma_noise,var0,"o")
title("Var (m0)")

nexttile
plot(sigma_noise,moment1,"o");
title("Mean (m1)")

nexttile
plot(sigma_noise,var1,"o")
title("Var (m1)")

nexttile
plot(sigma_noise,moment2,"o");
title("Mean (m2)")

nexttile
plot(sigma_noise,var2,"o")
title("Var (m2)")

nexttile
plot(sigma_noise,moment3,"o");
title("Mean (m3)")

nexttile
plot(sigma_noise,var3,"o")
title("Var (m3)")

%% Varying lambda
clear all

tic
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = linspace(10e6,20e9,50); % randomly chosen arrival rate lambda 10e9 arrivals per second
n = 30;
sigma_noise = sqrt(0.28e-9);
order = 3;
for k = 1:length(lambda)
    [P Pv alpha tau] = turin_sim_alpha_cartesian_form_claus(lambda(k),G0,T,n,sigma_noise);
    for l = 1:n
        moment0bef(l) = abs(alpha(l,:)).^2*tau(:,l).^0;
        moment0(k) = mean(moment0bef);
        var0(k) = var(moment0bef);
        
        moment1bef(l) = abs(alpha(l,:)).^2*tau(:,l).^1;
        moment1(k) = mean(moment1bef);
        var1(k) = var(moment1bef);
        
        moment2bef(l) = abs(alpha(l,:)).^2*tau(:,l).^2;
        moment2(k) = mean(moment2bef);
        var2(k) = var(moment2bef);
        
        moment3bef(l) = abs(alpha(l,:)).^2*tau(:,l).^3;
        moment3(k) = mean(moment3bef);
        var3(k) = var(moment3bef);
    end
    k
end    
toc
figure(5)
plotT = tiledlayout(4,2)
title(plotT,"Summary statistics varying lambda")

nexttile
plot(lambda,moment0,"o")
title("Mean (m0)")

nexttile
plot(lambda,var0,"o")
title("Var (m0)")

nexttile
plot(lambda,moment1,"o");
title("Mean (m1)")

nexttile
plot(lambda,var1,"o")
title("Var (m1)")

nexttile
plot(lambda,moment2,"o");
title("Mean (m2)")

nexttile
plot(lambda,var2,"o")
title("Var (m2)")

nexttile
plot(lambda,moment3,"o");
title("Mean (m3)")

nexttile
plot(lambda,var3,"o")
title("Var (m3)")
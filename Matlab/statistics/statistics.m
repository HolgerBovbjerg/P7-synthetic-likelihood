% Matlab script for finding summary statistics from samples from the Turin
% model. In eact section one parameter is varied while the others are held
% fixed. 50 different values of each parameter is simulated

% The temporal moments are found using equation 3 in the paper "Parameter
% estimation for stochastic channel models using temporal moments"

% When varying a value, n iterations are being done. The mean and variance
% of these iterations are then found as our summary statistics.


%% Varying T
clear all

%T = linspace(7.8e-11,7.8e-9,50); % Reverberation time, 
a = 7.8e-11; b = 7.8e-7,Nl = 50;
Tvary = a + (b-a).*rand(Nl,1);
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
n = 20; % Number of iterations 

sigma_noise = sqrt(0.28e-9); % noise variance
for k = 1:length(Tvary)
    [P Pv alpha tau, t] = turin_sim_alpha_cartesian_form_claus(lambda,G0,Tvary(k),n,sigma_noise);

    for i = 1:n
        temporal0T(k,i) = trapz(t,(t.^0.*Pv(:,i)'));
        temporal1T(k,i) = trapz(t,(t.^1.*Pv(:,i)'));
        temporal2T(k,i) = trapz(t,(t.^2.*Pv(:,i)'));
        temporal3T(k,i) = trapz(t,(t.^3.*Pv(:,i)'));
    end
    k
end

 
 figure(2)
plotmean =  tiledlayout(4,1)
title(plotmean,"Mean of the temporal moments while varying T") 
nexttile
plot(Tvary,mean(temporal0T,2),"o")
title("Mean (m0)")

 nexttile
plot(Tvary,mean(temporal1T,2),"o") 
title("Mean (m1)")

 nexttile
plot(Tvary,mean(temporal2T,2),"o") 
title("Mean (m2)")

 nexttile
plot(Tvary,mean(temporal3T,2),"o")
title("Mean (m3)")

figure(3)
plotvar = tiledlayout(4,1)
title(plotvar,"Variance of the temporal moments while varying T")

nexttile
plot(Tvary,var(temporal0T'),"o")
title("Var (m0)")

nexttile
plot(Tvary,var(temporal1T'),"o")
title("Var (m1)")

nexttile
plot(Tvary,var(temporal2T'),"o")
title("Var (m2)")

nexttile
plot(Tvary,var(temporal3T'),"o")
title("Var (m3)")


%% Varying G0
%clear all

%T = linspace(7.8e-11,7.8e-9,50); % Reverberation time, 
T = 7.8e-9;
a = db2pow(-110); b = db2pow(-60),Nl = 50;
G0vary = a + (b-a).*rand(Nl,1);
% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
n = 20; % Number of iterations 

sigma_noise = sqrt(0.28e-9); % noise variance
for k = 1:length(G0vary)
    [P Pv alpha tau, t] = turin_sim_alpha_cartesian_form_claus(lambda,G0vary(k),T,n,sigma_noise);

    for i = 1:n
        temporal0G0(k,i) = trapz(t,(t.^0.*Pv(:,i)'));
        temporal1G0(k,i) = trapz(t,(t.^1.*Pv(:,i)'));
        temporal2G0(k,i) = trapz(t,(t.^2.*Pv(:,i)'));
        temporal3G0(k,i) = trapz(t,(t.^3.*Pv(:,i)'));
    end
    k
end

 
 figure(4)
plotmean =  tiledlayout(4,1)
title(plotmean,"Mean of the temporal moments while varying G0") 
nexttile
plot(G0vary,mean(temporal0G0,2),"o")
title("Mean (m0)")

 nexttile
plot(G0vary,mean(temporal1G0,2),"o") 
title("Mean (m1)")

 nexttile
plot(G0vary,mean(temporal2G0,2),"o") 
title("Mean (m2)")

 nexttile
plot(G0vary,mean(temporal3G0,2),"o")
title("Mean (m3)")

figure(5)
plotvar = tiledlayout(4,1)
title(plotvar,"Variance of the temporal moments while varying G0")

nexttile
plot(G0vary,var(temporal0G0'),"o")
title("Var (m0)")

nexttile
plot(G0vary,var(temporal1G0'),"o")
title("Var (m1)")

nexttile
plot(G0vary,var(temporal2G0'),"o")
title("Var (m2)")

nexttile
plot(G0vary,var(temporal3G0'),"o")
title("Var (m3)")




%% Varying sigma_noise
%clear all

%T = linspace(7.8e-11,7.8e-9,50); % Reverberation time, 
T = 7.8e-9;
G0 = db2pow(-83.9);
a = sqrt(5e-10); b = sqrt(5e-9),Nl = 50;
sigma_noise_vary = a + (b-a).*rand(Nl,1);
% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
n = 20; % Number of iterations 

for k = 1:length(sigma_noise_vary)
    [P Pv alpha tau, t] = turin_sim_alpha_cartesian_form_claus(lambda,G0,T,n,sigma_noise_vary(k));

    for i = 1:n
        temporal0sigma_noise(k,i) = trapz(t,(t.^0.*Pv(:,i)'));
        temporal1sigma_noise(k,i) = trapz(t,(t.^1.*Pv(:,i)'));
        temporal2sigma_noise(k,i) = trapz(t,(t.^2.*Pv(:,i)'));
        temporal3sigma_noise(k,i) = trapz(t,(t.^3.*Pv(:,i)'));
    end
    k
end

 
 figure(6)
plotmean =  tiledlayout(4,1)
title(plotmean,"Mean of the temporal moments while varying sigma_noise") 
nexttile
plot(sigma_noise_vary,mean(temporal0sigma_noise,2),"o")
title("Mean (m0)")

 nexttile
plot(sigma_noise_vary,mean(temporal1sigma_noise,2),"o") 
title("Mean (m1)")

 nexttile
plot(sigma_noise_vary,mean(temporal2sigma_noise,2),"o") 
title("Mean (m2)")

 nexttile
plot(sigma_noise_vary,mean(temporal3sigma_noise,2),"o")
title("Mean (m3)")

figure(7)
plotvar = tiledlayout(4,1)
title(plotvar,"Variance of the temporal moments while varying sigma_noise")

nexttile
plot(sigma_noise_vary,var(temporal0sigma_noise'),"o")
title("Var (m0)")

nexttile
plot(sigma_noise_vary,var(temporal1sigma_noise'),"o")
title("Var (m1)")

nexttile
plot(sigma_noise_vary,var(temporal2sigma_noise'),"o")
title("Var (m2)")

nexttile
plot(sigma_noise_vary,var(temporal3sigma_noise'),"o")
title("Var (m3)")




%% Varying lambda
%clear all

%T = linspace(7.8e-11,7.8e-9,50); % Reverberation time, 
T = 7.8e-9;
G0 = db2pow(-83.9);
a = 50e6; b = 50e7,Nl = 50;
lambdavary = a + (b-a).*rand(Nl,1)
% Time delay tau is a possion arrival process with mean delay lambda
n = 20; % Number of iterations 

sigma_noise = sqrt(0.28e-9); % noise variance
for k = 1:length(lambdavary)
    [P Pv alpha tau, t] = turin_sim_alpha_cartesian_form_claus(lambdavary(k),G0,T,n,sigma_noise);

    for i = 1:n
        temporal0lambda(k,i) = trapz(t,(t.^0.*Pv(:,i)'));
        temporal1lambda(k,i) = trapz(t,(t.^1.*Pv(:,i)'));
        temporal2lambda(k,i) = trapz(t,(t.^2.*Pv(:,i)'));
        temporal3lambda(k,i) = trapz(t,(t.^3.*Pv(:,i)'));
    end
    k
end

 
 figure(8)
plotmean =  tiledlayout(4,1)
title(plotmean,"Mean of the temporal moments while varying lambda") 
nexttile
plot(lambdavary,mean(temporal0lambda,2),"o")
title("Mean (m0)")

 nexttile
plot(lambdavary,mean(temporal1lambda,2),"o") 
title("Mean (m1)")

 nexttile
plot(lambdavary,mean(temporal2lambda,2),"o") 
title("Mean (m2)")

 nexttile
plot(lambdavary,mean(temporal3lambda,2),"o")
title("Mean (m3)")

figure(9)
plotvar = tiledlayout(4,1)
title(plotvar,"Variance of the temporal moments while varying lambda")

nexttile
plot(lambdavary,var(temporal0lambda'),"o")
title("Var (m0)")

nexttile
plot(lambdavary,var(temporal1lambda'),"o")
title("Var (m1)")

nexttile
plot(lambdavary,var(temporal2lambda'),"o")
title("Var (m2)")

nexttile
plot(lambdavary,var(temporal3lambda'),"o")
title("Var (m3)")













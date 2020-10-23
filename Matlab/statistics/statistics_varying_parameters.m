% Matlab script for finding summary statistics from samples from the Turin
% model. In eact section one parameter is varied while the others are held
% fixed. 50 different values of each parameter is simulated

% The temporal moments are found using equation 3 in the paper "Parameter
% estimation for stochastic channel models using temporal moments"

% When varying a value, n iterations are being done. The mean and variance
% of these iterations are then found as our summary statistics.

%% Varying T
%clear all

%T = linspace(7.8e-11,7.8e-9,50); % Reverberation time, 
a = 7.8e-11; b = 7.8e-7,Nl = 100;
Tvary = a + (b-a).*rand(Nl,1);
B = 4e9; % Bandwidth of signal: 4 GHz
n = 200; % Number of iterations 

statsmatrixT = zeros(Nl,8);


G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second



sigma_N = sqrt(0.28e-9); % noise variance

%temporal0T = zeros(length(Tvary),n);
% for k = 1:length(Tvary)
%         [Pv, t] = sim_turin_matrix_gpu(n, B, 801, Tvary(k), G0, lambda, sigma_N);
% 
%     for i = 1:n
%         temporal0T(k,i) = trapz(t,(t.^0.*Pv(:,i)));
%         temporal1T(k,i) = trapz(t,(t.^1.*Pv(:,i)));
%         temporal2T(k,i) = trapz(t,(t.^2.*Pv(:,i)));
%         temporal3T(k,i) = trapz(t,(t.^3.*Pv(:,i)));
%     end
%     k
% end


 figure(2)
plotmean =  tiledlayout(4,1)
title(plotmean,"Mean of the temporal moments while varying T") 
nexttile
statsmatrixT(:,1) = mean(temporal0T,2);
plot(Tvary,statsmatrixT(:,1),"o")
title("Mean (m0)")

 nexttile
statsmatrixT(:,2) = mean(temporal1T,2);
plot(Tvary,statsmatrixT(:,2),"o")
title("Mean (m1)")

 nexttile
statsmatrixT(:,3) = mean(temporal2T,2);
plot(Tvary,statsmatrixT(:,3),"o")
title("Mean (m2)")

 nexttile
statsmatrixT(:,4) = mean(temporal3T,2);
plot(Tvary,statsmatrixT(:,4),"o")
title("Mean (m3)")


figure(3)
plotvar = tiledlayout(4,1)
title(plotvar,"Variance of the temporal moments while varying T")

nexttile
statsmatrixT(:,5) = var(temporal0T');
plot(Tvary,statsmatrixT(:,5),"o")
title("Var (m0)")

nexttile
statsmatrixT(:,6) = var(temporal1T');
plot(Tvary,statsmatrixT(:,6),"o")
title("Var (m1)")

nexttile
statsmatrixT(:,7) = var(temporal2T');
plot(Tvary,statsmatrixT(:,7),"o")
title("Var (m2)")

nexttile
statsmatrixT(:,8) = var(temporal3T');
plot(Tvary,statsmatrixT(:,8),"o")
title("Var (m3)")


%% Varying G0
%clear all

%T = linspace(7.8e-11,7.8e-9,50); % Reverberation time, 
T = 7.8e-9;
a = db2pow(-110); b = db2pow(-60),Nl = 100;
G0vary = a + (b-a).*rand(Nl,1);
% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
n = 100; % Number of iterations 
B = 4e9; % Bandwidth of signal: 4 GHz
statsmatrixG0 = zeros(Nl,8);

sigma_N = sqrt(0.28e-9); % noise variance
for k = 1:length(G0vary)
        [Pv, t] = sim_turin_matrix_gpu(n, B, 801, T, G0vary(k), lambda, sigma_N);

    for i = 1:n
        temporal0G0(k,i) = trapz(t,(t.^0.*Pv(:,i)));
        temporal1G0(k,i) = trapz(t,(t.^1.*Pv(:,i)));
        temporal2G0(k,i) = trapz(t,(t.^2.*Pv(:,i)));
        temporal3G0(k,i) = trapz(t,(t.^3.*Pv(:,i)));
    end
    k
end

 
%  figure(4)
% plotmean =  tiledlayout(4,1)
% title(plotmean,"Mean of the temporal moments while varying G0") 
% nexttile
% plot(G0vary,mean(temporal0G0,2),"o")
% title("Mean (m0)")
% 
%  nexttile
% plot(G0vary,mean(temporal1G0,2),"o") 
% title("Mean (m1)")
% 
%  nexttile
% plot(G0vary,mean(temporal2G0,2),"o") 
% title("Mean (m2)")
% 
%  nexttile
% plot(G0vary,mean(temporal3G0,2),"o")
% title("Mean (m3)")
% 
% figure(5)
% plotvar = tiledlayout(4,1)
% title(plotvar,"Variance of the temporal moments while varying G0")
% 
% nexttile
% plot(G0vary,var(temporal0G0'),"o")
% title("Var (m0)")
% 
% nexttile
% plot(G0vary,var(temporal1G0'),"o")
% title("Var (m1)")
% 
% nexttile
% plot(G0vary,var(temporal2G0'),"o")
% title("Var (m2)")
% 
% nexttile
% plot(G0vary,var(temporal3G0'),"o")
% title("Var (m3)")
% 

figure(8)
plotmean =  tiledlayout(4,1)
title(plotmean,"Mean of the temporal moments while varying G0") 
nexttile
statsmatrixG0(:,1) = mean(temporal0G0,2);
plot(G0vary,statsmatrixG0(:,1),"o")
title("Mean (m0)")

 nexttile
statsmatrixG0(:,2) = mean(temporal1G0,2);
plot(G0vary,statsmatrixG0(:,2),"o")
title("Mean (m1)")

 nexttile
statsmatrixG0(:,3) = mean(temporal2G0,2);
plot(G0vary,statsmatrixG0(:,3),"o")
title("Mean (m2)")

 nexttile
statsmatrixG0(:,4) = mean(temporal3G0,2);
plot(G0vary,statsmatrixG0(:,4),"o")
title("Mean (m3)")


figure(9)
plotvar = tiledlayout(4,1)
title(plotvar,"Variance of the temporal moments while varying G0")

nexttile
statsmatrixG0(:,5) = var(temporal0G0');
plot(G0vary,statsmatrixG0(:,5),"o")
title("Var (m0)")

nexttile
statsmatrixG0(:,6) = var(temporal1G0');
plot(G0vary,statsmatrixG0(:,6),"o")
title("Var (m1)")

nexttile
statsmatrixG0(:,7) = var(temporal2G0');
plot(G0vary,statsmatrixG0(:,7),"o")
title("Var (m2)")

nexttile
statsmatrixG0(:,8) = var(temporal3G0');
plot(G0vary,statsmatrixG0(:,8),"o")
title("Var (m3)")


%% Varying sigma_noise
%clear all

%T = linspace(7.8e-11,7.8e-9,50); % Reverberation time, 
T = 7.8e-9;
G0 = db2pow(-83.9);
a = sqrt(5e-10); b = sqrt(5e-9),Nl = 100;
sigma_noise_vary = a + (b-a).*rand(Nl,1);
% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
n = 100; % Number of iterations 
B = 4e9; % Bandwidth of signal: 4 GHz
statsmatrixsigma = zeros(Nl,8);
for k = 1:length(sigma_noise_vary)
        [Pv, t] = sim_turin_matrix_gpu(n, B, 801, T, G0, lambda, sigma_noise_vary(k));

    for i = 1:n
        temporal0sigma_noise(k,i) = trapz(t,(t.^0.*Pv(:,i)));
        temporal1sigma_noise(k,i) = trapz(t,(t.^1.*Pv(:,i)));
        temporal2sigma_noise(k,i) = trapz(t,(t.^2.*Pv(:,i)));
        temporal3sigma_noise(k,i) = trapz(t,(t.^3.*Pv(:,i)));
    end
    k
end

%  
%  figure(6)
% plotmean =  tiledlayout(4,1)
% title(plotmean,"Mean of the temporal moments while varying sigma_noise") 
% nexttile
% plot(sigma_noise_vary,mean(temporal0sigma_noise,2),"o")
% title("Mean (m0)")
% 
%  nexttile
% plot(sigma_noise_vary,mean(temporal1sigma_noise,2),"o") 
% title("Mean (m1)")
% 
%  nexttile
% plot(sigma_noise_vary,mean(temporal2sigma_noise,2),"o") 
% title("Mean (m2)")
% 
%  nexttile
% plot(sigma_noise_vary,mean(temporal3sigma_noise,2),"o")
% title("Mean (m3)")
% 
% figure(7)
% plotvar = tiledlayout(4,1)
% title(plotvar,"Variance of the temporal moments while varying sigma_noise")
% 
% nexttile
% plot(sigma_noise_vary,var(temporal0sigma_noise'),"o")
% title("Var (m0)")
% 
% nexttile
% plot(sigma_noise_vary,var(temporal1sigma_noise'),"o")
% title("Var (m1)")
% 
% nexttile
% plot(sigma_noise_vary,var(temporal2sigma_noise'),"o")
% title("Var (m2)")
% 
% nexttile
% plot(sigma_noise_vary,var(temporal3sigma_noise'),"o")
% title("Var (m3)")
% 

figure(8)
plotmean =  tiledlayout(4,1)
title(plotmean,"Mean of the temporal moments while varying \sigma") 
nexttile
statsmatrixsigma(:,1) = mean(temporal0sigma_noise,2);
plot(sigma_noise_vary,statsmatrixsigma(:,1),"o")
title("Mean (m0)")

 nexttile
statsmatrixsigma(:,2) = mean(temporal1sigma_noise,2);
plot(sigma_noise_vary,statsmatrixsigma(:,2),"o")
title("Mean (m1)")

 nexttile
statsmatrixsigma(:,3) = mean(temporal2sigma_noise,2);
plot(sigma_noise_vary,statsmatrixsigma(:,3),"o")
title("Mean (m2)")

 nexttile
statsmatrixsigma(:,4) = mean(temporal3sigma_noise,2);
plot(sigma_noise_vary,statsmatrixsigma(:,4),"o")
title("Mean (m3)")


figure(9)
plotvar = tiledlayout(4,1)
title(plotvar,"Variance of the temporal moments while varying \sigma")

nexttile
statsmatrixsigma(:,5) = var(temporal0sigma_noise');
plot(sigma_noise_vary,statsmatrixsigma(:,5),"o")
title("Var (m0)")

nexttile
statsmatrixsigma(:,6) = var(temporal1sigma_noise');
plot(sigma_noise_vary,statsmatrixsigma(:,6),"o")
title("Var (m1)")

nexttile
statsmatrixsigma(:,7) = var(temporal2sigma_noise');
plot(sigma_noise_vary,statsmatrixsigma(:,7),"o")
title("Var (m2)")

nexttile
statsmatrixsigma(:,8) = var(temporal3sigma_noise');
plot(sigma_noise_vary,statsmatrixsigma(:,8),"o")
title("Var (m3)")


%% Varying lambda
%clear all

%T = linspace(7.8e-11,7.8e-9,50); % Reverberation time, 
T = 7.8e-9;
G0 = db2pow(-83.9);
a = 50e6; b = 50e7,Nl = 100;
lambdavary = a + (b-a).*rand(Nl,1)
% Time delay tau is a possion arrival process with mean delay lambda
n = 100; % Number of iterations 
B = 4e9; % Bandwidth of signal: 4 GHz
statsmatrixLambda = zeros(Nl,8);

sigma_N = sqrt(0.28e-9); % noise variance
for k = 1:length(lambdavary)
        [Pv, t] = sim_turin_matrix_gpu(n, B, 801, T, G0, lambdavary(k), sigma_N);

    for i = 1:n
        temporal0lambda(k,i) = trapz(t,(t.^0.*Pv(:,i)));
        temporal1lambda(k,i) = trapz(t,(t.^1.*Pv(:,i)));
        temporal2lambda(k,i) = trapz(t,(t.^2.*Pv(:,i)));
        temporal3lambda(k,i) = trapz(t,(t.^3.*Pv(:,i)));
    end
    k
end

 
%  figure(8)
% plotmean =  tiledlayout(4,1)
% title(plotmean,"Mean of the temporal moments while varying lambda") 
% nexttile
% plot(lambdavary,mean(temporal0lambda,2),"o")
% title("Mean (m0)")
% 
%  nexttile
% plot(lambdavary,mean(temporal1lambda,2),"o") 
% title("Mean (m1)")
% 
%  nexttile
% plot(lambdavary,mean(temporal2lambda,2),"o") 
% title("Mean (m2)")
% 
%  nexttile
% plot(lambdavary,mean(temporal3lambda,2),"o")
% title("Mean (m3)")
% 
% figure(9)
% plotvar = tiledlayout(4,1)
% title(plotvar,"Variance of the temporal moments while varying lambda")
% 
% nexttile
% plot(lambdavary,var(temporal0lambda'),"o")
% title("Var (m0)")
% 
% nexttile
% plot(lambdavary,var(temporal1lambda'),"o")
% title("Var (m1)")
% 
% nexttile
% plot(lambdavary,var(temporal2lambda'),"o")
% title("Var (m2)")
% 
% nexttile
% plot(lambdavary,var(temporal3lambda'),"o")
% title("Var (m3)")
%%
 figure(8)
plotmean =  tiledlayout(4,1)
title(plotmean,"Mean of the temporal moments while varying lambda") 
nexttile
statsmatrixLambda(:,1) = mean(temporal0lambda,2);
plot(lambdavary,statsmatrixLambda(:,1),"o")
title("Mean (m0)")

 nexttile
statsmatrixLambda(:,2) = mean(temporal1lambda,2);
plot(lambdavary,statsmatrixLambda(:,2),"o")
title("Mean (m1)")

 nexttile
statsmatrixLambda(:,3) = mean(temporal2lambda,2);
plot(lambdavary,statsmatrixLambda(:,3),"o")
title("Mean (m2)")

 nexttile
statsmatrixLambda(:,4) = mean(temporal3lambda,2);
plot(lambdavary,statsmatrixLambda(:,4),"o")
title("Mean (m3)")


figure(9)
plotvar = tiledlayout(4,1)
title(plotvar,"Variance of the temporal moments while varying T")

nexttile
statsmatrixLambda(:,5) = var(temporal0lambda');
plot(lambdavary,statsmatrixLambda(:,5),"o")
title("Var (m0)")

nexttile
statsmatrixLambda(:,6) = var(temporal1lambda');
plot(lambdavary,statsmatrixLambda(:,6),"o")
title("Var (m1)")

nexttile
statsmatrixLambda(:,7) = var(temporal2lambda');
plot(lambdavary,statsmatrixLambda(:,7),"o")
title("Var (m2)")

nexttile
statsmatrixLambda(:,8) = var(temporal3lambda');
plot(lambdavary,statsmatrixLambda(:,8),"o")
title("Var (m3)")


%%

T = 7.8e-9;
G0 = db2pow(-83.9);

a = 7.8e-11; b = 7.8e-7,Nl = 100;
Tvary = a + (b-a).*rand(Nl,1);

a = 50e6; b = 50e7,Nl = 100;
lambdavary = a + (b-a).*rand(Nl,1)


a = sqrt(5e-10); b = sqrt(5e-9),Nl = 100;
sigma_noise_vary = a + (b-a).*rand(Nl,1);


a = db2pow(-110); b = db2pow(-60),Nl = 100;
G0vary = a + (b-a).*rand(Nl,1);



% Time delay tau is a possion arrival process with mean delay lambda
n = 100; % Number of iterations 
B = 4e9; % Bandwidth of signal: 4 GHz
statsmatrixeverything = zeros(Nl,8);

%sigma_N = sqrt(0.28e-9); % noise variance
for k = 1:length(lambdavary)
        [Pv, t] = sim_turin_matrix_gpu(n, B, 801, Tvary(k), G0vary(k), lambdavary(k), sigma_noise_vary(k));

    for i = 1:n
        temporal0e(k,i) = trapz(t,(t.^0.*Pv(:,i)));
        temporal1e(k,i) = trapz(t,(t.^1.*Pv(:,i)));
        temporal2e(k,i) = trapz(t,(t.^2.*Pv(:,i)));
        temporal3e(k,i) = trapz(t,(t.^3.*Pv(:,i)));
    end
    k
end


 figure(10)
plotmean =  tiledlayout(4,1)
title(plotmean,"Mean of the temporal moments while varying everything") 
nexttile
statsmatrixeverything(:,1) = mean(temporal0e,2);
plot(lambdavary,statsmatrixeverything(:,1),"o")
title("Mean (m0)")

 nexttile
statsmatrixeverything(:,2) = mean(temporal1e,2);
plot(lambdavary,statsmatrixeverything(:,2),"o")
title("Mean (m1)")

 nexttile
statsmatrixeverything(:,3) = mean(temporal2e,2);
plot(lambdavary,statsmatrixeverything(:,3),"o")
title("Mean (m2)")

 nexttile
statsmatrixeverything(:,4) = mean(temporal3e,2);
plot(lambdavary,statsmatrixeverything(:,4),"o")
title("Mean (m3)")


figure(9)
plotvar = tiledlayout(4,1)
title(plotvar,"Variance of the temporal moments while varying T")

nexttile
statsmatrixeverything(:,5) = var(temporal0e');
plot(lambdavary,statsmatrixeverything(:,5),"o")
title("Var (m0)")

nexttile
statsmatrixeverything(:,6) = var(temporal1e');
plot(lambdavary,statsmatrixeverything(:,6),"o")
title("Var (m1)")

nexttile
statsmatrixeverything(:,7) = var(temporal2e');
plot(lambdavary,statsmatrixeverything(:,7),"o")
title("Var (m2)")

nexttile
statsmatrixeverything(:,8) = var(temporal3e');
plot(lambdavary,statsmatrixeverything(:,8),"o")
title("Var (m3)")


%% Plot correlation map

Tcoeffs = corrcoef(statsmatrixT);
lambdacoeffs = corrcoef(statsmatrixLambda);
sigmacoeffs = corrcoef(statsmatrixsigma);
G0coeffs = corrcoef(statsmatrixG0);
eveythingcoeffs = corrrcoef(statsmatrixeverything);

Clims = [0.5 1];
tiledlayout(5,2)

nexttile
imagesc(Tcoeffs, Clims);
colorbar
title("Correlation varying T")

nexttile
imagesc(Tcoeffs);
colorbar
title("Correlation varying T")

nexttile
imagesc(lambdacoeffs, Clims)
colorbar
title("Correlation varying \lambda")

nexttile
imagesc(lambdacoeffs)
colorbar
title("Correlation varying \lambda")

nexttile
imagesc(sigmacoeffs, Clims);
colorbar
title("Correlation varying \sigma")

nexttile
imagesc(sigmacoeffs);
colorbar
title("Correlation varying \sigma")

nexttile
imagesc(G0coeffs, Clims);
colorbar
title("Correlation varying G0")

nexttile
imagesc(G0coeffs);
colorbar
title("Correlation varying G0")

nexttile
imagesc(eveythingcoeffs, Clims);
colorbar
title("Correlation varying everything")

nexttile
imagesc(eveythingcoeffs);
colorbar
title("Correlation varying everything")

%%
eveythingcoeffs = corrcoef(statsmatrixeverything);

imagesc(eveythingcoeffs);
colorbar
title("Correlation varying everything")
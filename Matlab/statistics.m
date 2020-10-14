% %%
% 
% clear all
% tic
% T = 7.8e-9;
% %T = linspace(5e-9, 10e-9, 50); % Reverberation time: 7.8 ns
% G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
% 
% % Time delay tau is a possion arrival process with mean delay lambda
% %lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second
% lambda = linspace(5e7, 15e10, 75);
% n = 3;
% % mean0 = zeros(1,50);
% % for i = 1:length(T)
% %     P = (turin_sim_alpha_cartesian_form_claus(lambda,G0,T(i),n));
% %     mean0(i) = mean(P);
% %     i
% % end
% mean0 = zeros(1,length(lambda));
% mean1 = zeros(1,length(lambda));
% mean2 = zeros(1,length(lambda));
% mean3 = zeros(1,length(lambda));
% for i = 1:length(lambda)
%     [P Pv] = (turin_sim_alpha_cartesian_form_claus(lambda(i),G0,T,n));
%     moment1 = moment(P,1);
%     moment2 = moment(P,2);
%     moment3 = moment(P,3);
%     mean0(i) = mean(P);
%     mean1(i) = mean(moment1);
%     mean2(i) = mean(moment2);
%     mean3(i) = mean(moment3);
%     
%     var0(i) = var(P);
%     var1(i) = var(moment(Pv,1));
%     var2(i) = var(moment(Pv,2));
%     var3(i) = var(moment(Pv,3));
%     i
% end
% toc
% figure(2)
% tiledlayout(4,2)
% nexttile
% plot(lambda,mean0,"o")
% title("mean(m0)")
% nexttile
% plot(lambda,var0,"o")
% title("var(m0)")
% nexttile
% plot(lambda,mean1,"o")
% title("mean(m1)")
% nexttile
% plot(lambda,var1,"o")
% title("var(m1)")
% nexttile
% plot(lambda,mean2,"o")
% title("mean(m2)")
% nexttile
% plot(lambda,var2,"o")
% title("var(m2)")
% nexttile
% plot(lambda,mean3,"o")
% title("mean(m3)")
% nexttile
% plot(lambda,var3,"o")
% title("var(m3)")
% 
% %%
% 
% clear all
% 
% tic
% T = 7.8e-9;
% %T = linspace(5e-9, 10e-9, 50); % Reverberation time: 7.8 ns
% G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
% 
% % Time delay tau is a possion arrival process with mean delay lambda
% %lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second
% lambda = linspace(5e7, 10e8, 50);
% n = 10;
% %moment1 = zeros(length(lambda))';
% 
% order = 3;
% for k = 1:length(lambda)
%     [P Pv alpha tau] = turin_sim_alpha_cartesian_form_claus(lambda(k),G0,T,n);
%     for l = 1:n
%         moment1bef(l) = abs(alpha(l,:)).^2*tau(:,l).^order;
%         moment1(k) = mean(moment1bef);
%     end
%     k
% end    
% toc
% figure(2)
% plot(lambda,moment1,"o");


%% Varying T
clear all

tic
T = linspace(7.8e-11,5e-8,50);
%T = linspace(5e-9, 10e-9, 50); % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
%lambda = linspace(5e7, 10e8, 50);
n = 15;
%moment1 = zeros(length(lambda))';
sigma_noise = sqrt(28e-9);
order = 3;
for k = 1:length(T)
    [P Pv alpha tau] = turin_sim_alpha_cartesian_form_claus(lambda,G0,T(k),n,sigma_noise);
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
%T = linspace(7.8e-11,5e-8,50);
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = linspace(db2pow(-110), db2pow(-60), 50); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
%lambda = linspace(5e7, 10e8, 50);
n = 15;
%moment1 = zeros(length(lambda))';
sigma_noise = sqrt(28e-9);
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
%T = linspace(7.8e-11,5e-8,50);
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
%lambda = linspace(5e7, 10e8, 50);
n = 15;
%moment1 = zeros(length(lambda))';
sigma_noise = linspace(sqrt(28e-11),sqrt(28e-7),50);
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
%T = linspace(7.8e-11,5e-8,50);
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = linspace(10e6,20e9,50); % randomly chosen arrival rate lambda 10e9 arrivals per second
%lambda = linspace(5e7, 10e8, 50);
n = 10;
%moment1 = zeros(length(lambda))';
sigma_noise = sqrt(28e-9);
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


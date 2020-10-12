%%

clear all
tic
T = 7.8e-9;
%T = linspace(5e-9, 10e-9, 50); % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
%lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second
lambda = linspace(5e7, 15e10, 75);
n = 15;
% mean0 = zeros(1,50);
% for i = 1:length(T)
%     P = (turin_sim_alpha_cartesian_form_claus(lambda,G0,T(i),n));
%     mean0(i) = mean(P);
%     i
% end
mean0 = zeros(1,length(lambda));
mean1 = zeros(1,length(lambda));
mean2 = zeros(1,length(lambda));
mean3 = zeros(1,length(lambda));
for i = 1:length(lambda)
    [P Pv] = (turin_sim_alpha_cartesian_form_claus(lambda(i),G0,T,n));
    moment1 = moment(P,1);
    moment2 = moment(P,2);
    moment3 = moment(P,3);
    mean0(i) = mean(P);
    mean1(i) = mean(moment1);
    mean2(i) = mean(moment2);
    mean3(i) = mean(moment3);
    
    var0(i) = var(P);
    var1(i) = var(moment(Pv,1));
    var2(i) = var(moment(Pv,2));
    var3(i) = var(moment(Pv,3));
    i
end
toc
figure(2)
tiledlayout(4,2)
nexttile
plot(lambda,mean0,"o")
title("mean(m0)")
nexttile
plot(lambda,var0,"o")
title("var(m0)")
nexttile
plot(lambda,mean1,"o")
title("mean(m1)")
nexttile
plot(lambda,var1,"o")
title("var(m1)")
nexttile
plot(lambda,mean2,"o")
title("mean(m2)")
nexttile
plot(lambda,var2,"o")
title("var(m2)")
nexttile
plot(lambda,mean3,"o")
title("mean(m3)")
nexttile
plot(lambda,var3,"o")
title("var(m3)")

%%

clear all

tic
T = 7.8e-9;
%T = linspace(5e-9, 10e-9, 50); % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
%lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second
lambda = linspace(5e7, 10e8, 50);
n = 10;
%moment1 = zeros(length(lambda))';

order = 3;
for k = 1:length(lambda)
    [P Pv alpha tau] = turin_sim_alpha_cartesian_form_claus(lambda(k),G0,T,n);
    for l = 1:n
        moment1bef(l) = abs(alpha(l,:)).^2*tau(:,l).^order;
        moment1(k) = mean(moment1bef);
    end
    k
end    
toc
figure(2)
plot(lambda,moment1,"o");


%%
clear all

tic
T = linspace(7.8e-11,5e-7,50);
%T = linspace(5e-9, 10e-9, 50); % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e8; % randomly chosen arrival rate lambda 10e9 arrivals per second
%lambda = linspace(5e7, 10e8, 50);
n = 10;
%moment1 = zeros(length(lambda))';

order = 3;
for k = 1:length(T)
    [P Pv alpha tau] = turin_sim_alpha_cartesian_form_claus(lambda,G0,T(k),n);
    for l = 1:n
        moment1bef(l) = abs(alpha(l,:)).^2*tau(:,l).^1;
        moment1(k) = mean(moment1bef);
        moment2bef(l) = abs(alpha(l,:)).^2*tau(:,l).^2;
        moment2(k) = mean(moment2bef);
        moment3bef(l) = abs(alpha(l,:)).^2*tau(:,l).^3;
        moment3(k) = mean(moment3bef);
    end
    k
end    
toc
figure(2)
tiledlayout(3,1)
nexttile
plot(T,moment1,"o");
nexttile
plot(T,moment2,"o");
nexttile
plot(T,moment3,"o");
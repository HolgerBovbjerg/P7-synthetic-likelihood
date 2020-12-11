%Test for increasing L

clear all

% Data are randomly drawn from the prior range of
% "Prior_data_small_prior_min_max_values.mat"
load("thetas_for_turin_tests.mat")

Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

sims = 100;
L = 30000;
N = 100;
PV_temp = zeros(801,sims);
s_sim = zeros(L,9);
s_mean = s_sim;
counter = 0;
for i = 1:1
    i
    parfor ii = 1:L
        counter = counter+1
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, turin_thetas_2);
        s_sim(ii,:) = create_statistics_new(Pv, t);
    end
    s_mean(i,:) = mean(s_sim);
end
figure(3)
plot(s_mean(:,1))

%plot(s_sim(:,1))

%%
figure(4)
tiledlayout(3,3)
set(0,'defaultTextInterpreter','latex');
bins=30;
% for p = 1:9
%    nexttile
%    hist(s_sim(:,p),bins)
% end

%---------------------------------------------------
nexttile
%hist(s_sim(:,1),bins)
histfit(s_sim(:,1),bins)
title('$mean(m_0)$')
set(gca,'XTick',[], 'YTick', [])

nexttile
%hist(s_sim(:,2),bins)
histfit(s_sim(:,2),bins)
title('$mean(m_1)$')
set(gca,'XTick',[], 'YTick', [])

nexttile
%hist(s_sim(:,3),bins)
histfit(s_sim(:,3),bins)
title('$mean(m_2)$')
set(gca,'XTick',[], 'YTick', [])

%---------------------------------------------------

nexttile
%hist(s_sim(:,4),bins)
histfit(s_sim(:,4),bins)
title('$Var(m_0)$')
set(gca,'XTick',[], 'YTick', [])

nexttile
%hist(s_sim(:,5),bins)
histfit(s_sim(:,5),bins)
title('$Var(m_1)$')
set(gca,'XTick',[], 'YTick', [])

nexttile
%hist(s_sim(:,6),bins)
histfit(s_sim(:,6),bins)
title('$Var(m_2)$')
set(gca,'XTick',[], 'YTick', [])

%---------------------------------------------------

nexttile
%hist(s_sim(:,7),bins)
histfit(s_sim(:,7),bins)
title('$Cov(m_0,m_1)$')
set(gca,'XTick',[], 'YTick', [])

nexttile
%hist(s_sim(:,8),bins)
histfit(s_sim(:,8),bins)
title('$Cov(m_0,m_2)$')
set(gca,'XTick',[], 'YTick', [])

nexttile
%hist(s_sim(:,9),bins)
histfit(s_sim(:,9),bins)
title('$Cov(m_1,m_2)$')
set(gca,'XTick',[], 'YTick', [])

%%
%Test for increasing N
clear all

% Data are randomly drawn from the prior range of
% "Prior_data_small_prior_min_max_values.mat"
load("thetas_for_turin_tests.mat")

Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

sims = 2000;
L = 1;
N = 2000;
PV_temp = zeros(801,sims);
s_sim = zeros(N,9);
s_mean = s_sim;
for i = 1:L
    for ii = 100:N
        ii
        [Pv, t] = sim_turin_matrix_gpu(ii, B, Ns, turin_thetas_2);
        s_sim(ii,:) = create_statistics_new(Pv, t);
    end
    s_mean(i,:) = mean(s_sim(1:i,:),1);
end

plot(s_mean(:,1))

plot(s_sim(100:N,1))

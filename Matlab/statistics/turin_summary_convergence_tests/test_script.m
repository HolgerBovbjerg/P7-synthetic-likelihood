%Test for increasing L

clear all

% Data are randomly drawn from the prior range of
% "Prior_data_small_prior_min_max_values.mat"
load("thetas_for_turin_tests.mat")

Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

sims = 2000;
L = 800;
N = 300;
PV_temp = zeros(801,sims);
s_sim = zeros(L,9);
s_mean = s_sim;
for i = 1:L
    i
    parfor ii = 1:i
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, turin_thetas_2);
        s_sim(ii,:) = create_statistics_new(Pv, t);
    end
    s_mean(i,:) = mean(s_sim(1:i,:),1);
end
figure(3)
plot(s_mean(:,1))

%plot(s_sim(:,1))


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

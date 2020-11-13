clear

load('Prior_data_medium_prior_min_max_values.mat')
load('Theta_true_values.mat')


N = 200; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

tiledlayout(9,4)

j = 1;
steps = 100;

% T
S = zeros(steps,9);
T = linspace(1e-9,20e-9);
G0 = theta_true(2);
lambda = theta_true(3);
sigma_N = theta_true(4);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T(step) G0 lambda sigma_N]);
    S(step,:) = create_statistics_new(Pv, t);    
end
ylabels =['mean(m_0)','mean(m_1)','mean(m_2)','var(m_0)','var(m_1)','var(m_2)','cov(m_0,m_1)','cov(m_0,m_2)','cov(m_1,m_2)'];
tiles = [1 5 9 13 17 21 25 29 33];
for i = 1:9
    nexttile(tiles(i))
    scatter(T,S(:,i))
    if i == 1
        title('T')
    end
    ylabel(ylabels(i))
end

% G0
S = zeros(steps,9);
T = theta_true(1);
G0db = linspace(-90, -70, steps);
lambda = theta_true(3);
sigma_N = theta_true(4);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T db2pow(G0db(step)) lambda sigma_N]);
    S(step,:) = create_statistics_new(Pv, t);    
end

for i = 1:9
    nexttile(tiles(i)+1)
    scatter(G0db,S(:,i))
    if i == 1
        title('G_0')
    end
end

% lambda
S = zeros(steps,9);
T = theta_true(1);
G0 = theta_true(2);
lambda = linspace(1e8,1e10);
sigma_N = theta_true(4);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T G0 lambda(step) sigma_N]);
    S(step,:) = create_statistics_new(Pv, t);    
end

for i = 1:9
    nexttile(tiles(i)+2)
    scatter(lambda,S(:,i))
    if i == 1
        title('\lambda')
    end
end

% sigma_N
S = zeros(steps,9);
T = theta_true(1);
G0 = theta_true(2);
lambda = theta_true(3);
sigma_N = linspace(sqrt(0.28e-10), sqrt(0.28e-8),steps);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [G0 T lambda sigma_N(step)]);
    S(step,:) = create_statistics_new(Pv, t);    
end

for i = 1:9
    nexttile(tiles(i)+3)
    scatter(sigma_N,S(:,i))
    if i == 1
        title('\sigma_N')
    end
end
% nexttile
% plot(S(:,2))
% nexttile
% plot(S(:,3))
% nexttile
% plot(S(:,4))
% nexttile
% plot(S(:,5))
% nexttile
% plot(S(:,6))
% nexttile
% plot(S(:,7))
% nexttile
% plot(S(:,8))
% nexttile
% plot(S(:,9))

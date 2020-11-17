clear

load('Prior_data_large_prior_min_max_values.mat')
load('Theta_true_values.mat')


N = 200; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

tiledlayout(9,4)

j = 1;
steps = 100;

% T
S_T = zeros(steps,9);
T_vary = linspace(prior(1,1),prior(1,2));
G0 = theta_true(2);
lambda = theta_true(3);
sigma_N = theta_true(4);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T_vary(step) G0 lambda sigma_N]);
    S(step,:) = create_statistics(Pv, t);    
end
ylabels =["\mu(m_0)","\mu(m_1)","\mu(m_2)","var(m_0)","var(m_1)","var(m_2)","cov(m_0,m_1)","cov(m_0,m_2)","cov(m_1,m_2)"];
tiles = [1 5 9 13 17 21 25 29 33];
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])

% G0
S = zeros(steps,9);
T = theta_true(1);
G0db = linspace(-90, -70, steps);
lambda = theta_true(3);
sigma_N = theta_true(4);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T db2pow(G0db(step)) lambda sigma_N]);
    S(step,:) = create_statistics(Pv, t);    
end

% lambda
S = zeros(steps,9);
T = theta_true(1);
G0 = theta_true(2);
lambda = linspace(prior(3,1),prior(3,2));
sigma_N = theta_true(4);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T G0 lambda(step) sigma_N]);
    S(step,:) = create_statistics(Pv, t);    
end



% sigma_N
S = zeros(steps,9);
T = theta_true(1);
G0 = theta_true(2);
lambda = theta_true(3);
sigma_N_var = linspace(prior(4,1),prior(4,2),steps);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [G0 T lambda sigma_N_var(step)]);
    S_sigma(step,:) = create_statistics(Pv, t);    
end


%% Plots
for i = 1:9
    nexttile(tiles(i)) 
    scatter(T,S(:,i))
    if i == 1
        title('T')
    end
    ylabel(ylabels(i))
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end
for i = 1:9
    nexttile(tiles(i)+1)
    scatter(G0db,S(:,i))
    if i == 1
        title('G_0')
    end
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end
for i = 1:9
    nexttile(tiles(i)+2)
    scatter(lambda,S(:,i))
    if i == 1
        title('\lambda')
    end
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end
for i = 1:9
    nexttile(tiles(i)+3)
    scatter(sigma_N,S(:,i))
    if i == 1
        title('\sigma_N')
    end
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end




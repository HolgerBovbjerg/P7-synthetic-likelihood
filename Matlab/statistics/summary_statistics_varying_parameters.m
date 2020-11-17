clear

load('Prior_data_large_prior_min_max_values.mat')
load('Theta_true_values.mat')


N = 200; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

steps = 100;

% T
S_T = zeros(steps,9);
T_vary = linspace(prior(1,1),prior(1,2));
G0 = theta_true(2);
lambda = theta_true(3);
sigma_N = theta_true(4);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T_vary(step) G0 lambda sigma_N]);
    S_T(step,:) = create_statistics(Pv, t);    
end

% G0
S_G0 = zeros(steps,9);
T = theta_true(1);
G0db_vary = linspace(-90, -70, steps);
lambda = theta_true(3);
sigma_N = theta_true(4);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T db2pow(G0db_vary(step)) lambda sigma_N]);
    S_G0(step,:) = create_statistics(Pv, t);    
end

% lambda
S_lambda = zeros(steps,9);
T = theta_true(1);
G0 = theta_true(2);
lambda_vary = linspace(prior(3,1),prior(3,2));
sigma_N = theta_true(4);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T G0 lambda_vary(step) sigma_N]);
    S_lambda(step,:) = create_statistics(Pv, t);    
end

% sigma_N
S_sigma_N = zeros(steps,9);
T = theta_true(1);
G0 = theta_true(2);
lambda = theta_true(3);
sigma_N_vary = linspace(prior(4,1),prior(4,2),steps);
parfor step = 1:steps
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [G0 T lambda sigma_N_vary(step)]);
    S_sigma_N(step,:) = create_statistics(Pv, t);    
end


%% Plots
t = tiledlayout(9,4,'TileSpacing','Compact','Padding','Compact');
% t.Units = 'inches';
% t.OuterPosition = [0 0 3.5 7];
% ylabels =["\mu(m_0)","\mu(m_1)","\mu(m_2)","var(m_0)","var(m_1)","var(m_2)","cov(m_0,m_1)","cov(m_0,m_2)","cov(m_1,m_2)"];
ylabels =["s(1)","s(2)","s(3)","s(4)","s(5)","s(6)","s(7)","s(8)","s(9)"];

tiles(1,:) = [1 5 9 13 17 21 25 29 33];
set(t,'DefaultTextFontSize',8)
for i = 1:(9)
    nexttile(tiles(i)) 
    h = scatter(T_vary/1e-8,S_T(:,i));
    set(h,'SizeData',10)
    if i == 1
        subtitle('T \times 10^{-8}')
    end
    ylabel(ylabels(i))
    if (i ~= 9)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
    end
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
%     set(gca,'FontSize',8)   
end
for i = 1:9
    nexttile(tiles(i)+1)
    h = scatter(G0db_vary,S_G0(:,i));
    set(h,'SizeData',10)
    if i == 1
        subtitle('G_0 [dB]')
    end
    if (i ~= 9)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
    end
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
%     set(gca,'FontSize',8)
end
for i = 1:9
    nexttile(tiles(i)+2)
    h = scatter(lambda_vary/1e9,S_T(:,i));
    set(h,'SizeData',10)
    if i == 1
        subtitle('\lambda \times 10^{9}')
    end
    if (i ~= 9)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
    end
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
%     set(gca,'FontSize',8)
end
for i = 1:9
    nexttile(tiles(i)+3)
    h = scatter(sigma_N_vary.^2*1e9,S_sigma_N(:,i));
    set(h,'SizeData',10)
    if i == 1
        subtitle('\sigma_N^2 \times 10^{-9}')
    end
    if (i ~= 9)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
    end
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
%     set(gca,'FontSize',8)
end

%% Export
exportgraphics(t,'summary_information.pdf')




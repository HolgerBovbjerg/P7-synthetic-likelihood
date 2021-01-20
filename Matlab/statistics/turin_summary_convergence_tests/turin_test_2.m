clear all

% Data are randomly drawn from the prior range of
% "Prior_data_small_prior_min_max_values.mat"
load("thetas_for_turin_tests.mat")

Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
Q = 10;
sims = 1000;

PV_temp = zeros(801,sims);
s_obs = zeros(9,sims,Q,4);


thetas = [turin_thetas_1 turin_thetas_2 turin_thetas_3 turin_thetas_4]; 
for k = 1:4
    tic
    for j = 1:Q
        j
        [Pv, t] = sim_turin_matrix_gpu(1000, B, Ns, thetas(:,k));
        for i = 2:sims
            %     [Pv, t] = sim_turin_matrix_gpu(1, B, Ns, turin_thetas_1);
            s_obs(:,i,j,k) = create_statistics(Pv(:,1:i), t);
        end
    end
    toc
end
%%


colors =['b' 'r' 'g' 'm'];
tiledlayout(3,1)

nexttile
hold on
for k = 1:4

for j = 1:Q
    plot(s_obs(7,2:end-9,j,k),colors(k))
end
end
title('mean(m0)')
set(gca,'ytick',[])

nexttile
hold on
for k = 1:4
for j = 1:Q
    plot(s_obs(8,2:end-9,j,k), colors(k))
end
end
title('Var(m1)')
set(gca,'ytick',[])

nexttile
hold on
for k = 1:4
for j = 1:Q
    plot(s_obs(9,2:end-9,j,k),colors(k))
end
end
title('Cov(m0,m1)')
set(gca,'ytick',[])

%%
figure(2)
tiledlayout(4,1)

nexttile
hold on
for w=1:4
plot(thetas(1,w),'Color',colors(w),'Marker','o')
end

nexttile
hold on
for w=1:4
plot(thetas(2,w),'Color',colors(w),'Marker','o')
end

nexttile
hold on
for w=1:4
plot(thetas(3,w),'Color',colors(w),'Marker','o')
end

nexttile
hold on
for w=1:4
plot(thetas(4,w),'Color',colors(w),'Marker','o')
end



%%
load("test_using_turin_thetas_1.mat")
data1 = s_obs;

load("test_using_turin_thetas_2.mat")
data2 = s_obs;

load("test_using_turin_thetas_3.mat")
data3 = s_obs;

load("test_using_turin_thetas_4.mat")
data4 = s_obs;

ss_obs = [data1; data2; data3; data4];
%%
figure(5)
hold on
tiledlayout(4,1)
for p = 1:4
    nexttile
    hold on
    for q = 1:4
        
        plot(1:1000,ss_obs(4*(q-1)+p,1:1000),'Linewidth',2)
        if q == 1
            title('$\bar{m}_0$','Interpreter','Latex')
        elseif q == 2
            title('$\bar{m}_1$','Interpreter','Latex')
        elseif q == 3
            title('$Var(m_0)$','Interpreter','Latex')
        elseif q == 4
            title('$Var(m_1)$','Interpreter','Latex')
        end
    end
end
%%
[covariance theta_curr] = find_cov_prior(prior);
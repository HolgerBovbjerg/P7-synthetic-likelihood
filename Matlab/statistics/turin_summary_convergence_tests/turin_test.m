clear all

% Data are randomly drawn from the prior range of
% "Prior_data_small_prior_min_max_values.mat"
load("thetas_for_turin_tests.mat")

Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

sims = 2000;

PV_temp = zeros(801,sims);
s_obs = zeros(9,sims);

    [Pv, t] = sim_turin_matrix_gpu(1, B, Ns, turin_thetas_1);
    PV_temp(:,1) = Pv;
for i = 2:sims
    i
    [Pv, t] = sim_turin_matrix_gpu(1, B, Ns, turin_thetas_1);
    PV_temp(:,i) = Pv;
    s_obs(:,i) = create_statistics(PV_temp(:,1:i), t);
end
%%



tiledlayout(3,1)

nexttile
plot(s_obs(1,2:end-9))
   title('mean(m0)')
   set(gca,'ytick',[])
nexttile
plot(s_obs(4,2:end-9))
   title('Var(m1)')
   set(gca,'ytick',[])
nexttile
plot(s_obs(7,2:end-9))
   title('Cov(m0,m1)')
   set(gca,'ytick',[])


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
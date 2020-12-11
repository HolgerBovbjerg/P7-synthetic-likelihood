

load('Prior_data_small_prior_min_max_values.mat')
N = 100;
a = prior(1,1); b = prior(1,2); steps = 500;
T_vary = a + (b-a).*rand(steps,1);
a = prior(2,1); b = prior(2,2);
G0_vary = a + (b-a).*rand(steps,1);
a = prior(3,1); b = prior(3,2);
lambda_vary = a + (b-a).*rand(steps,1);
a = prior(4,1); b = prior(4,2);
sigma_noise_vary = a + (b-a).*rand(steps,1);

S_T = zeros(steps,9);
parfor step = 1:steps
    step
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, [T_vary(step) G0_vary(step) lambda_vary(step) sigma_noise_vary(step)]);
    S_T(step,:) = create_statistics_new(Pv, t);    
end

%%
co = corrcoef(S_T);

imagesc(co)
set(gca,'TickLabelInterpreter', 'latex');
ax = gca;
title(ax,'Correlation between summary statistics')
ax.YTick = 1:9;
ax.YTickLabel = {'$mean(m_0)$','$mean(m_1)$','$mean(m_2)$','$Var(m_0)$','$Var(m_1)$','$Var(m_2)$','$Cov(m_0,m_1)$','$Cov(m_0,m_2)$', '$Cov(m_1,m_2)$'};

ax.XTick = 1:9;
ax.XTickLabel = {'$mean(m_0)$','$mean(m_1)$','$mean(m_2)$','$Var(m_0)$','$Var(m_1)$','$Var(m_2)$','$Cov(m_0,m_1)$','$Cov(m_0,m_2)$', '$Cov(m_1,m_2)$'};

colorbar
xtickangle(ax,45)


clear all

N = 200; % Number of Turin simulations

T_true = 7.8e-9;
G0_true = db2pow(-83.9); % Reverberation gain converted from dB to power - 4.0738e-9
lambda_true = 4e8;
sigma_N_true = sqrt(0.28e-9); % noise variance

theta_true = [T_true G0_true lambda_true sigma_N_true];
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

sims = 1000;

PV_temp = zeros(801,sims);
s_obs = zeros(4,sims);
% "Observed data for testing"
for i = 1:sims
    i
    [Pv, t] = sim_turin_matrix_gpu(1, B, Ns, theta_true);
    PV_temp(:,i) = Pv;
    s_obs(:,i) = create_statistics(PV_temp(:,1:i), t);
end
%%
tiledlayout(4,1)
for q = 1:4
   nexttile
   plot(1:i,s_obs(q,:),'Linewidth',2)
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
%%
[covariance theta_curr] = find_cov_prior(prior);
clear
N = 300;
Bw = 4e9;
Ns = 801;

load('MMSE_est_BSL.mat')
[Pv_BSL, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_est_BSL);
load('MMSE_est_ABC.mat')
[Pv_ABC, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_est_ABC);
load("Theta_true_values.mat") 
[Pv_true,t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_true);

t = t*1e9;
tiles = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(t,pow2db(mean(Pv_true,2)),'black',t,pow2db(mean(Pv_BSL,2)),'green')
legend('True','BSL est.')
xtickformat('%d ns')
ylabel('[dB]')
box off
nexttile
plot(t,pow2db(mean(Pv_true,2)),'black',t,pow2db(mean(Pv_ABC,2)),'blue')
legend('True','ABC est.')
xtickformat('%d ns')
ylabel('[dB]')
box off
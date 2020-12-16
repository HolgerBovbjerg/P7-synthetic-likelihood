clear
N = 300;
Bw = 4e9;
Ns = 801;

load('MMSE_est_BSL.mat')
[Pv_BSL, ~] = sim_turin_matrix_gpu(N, Bw, Ns, theta_est_BSL);
load('MMSE_est_ABC_iter1.mat')
[Pv_ABC_iter1, ~] = sim_turin_matrix_gpu(N, Bw, Ns, theta_est_ABC);
load('MMSE_est_ABC_iter20.mat')
[Pv_ABC_iter20, ~] = sim_turin_matrix_gpu(N, Bw, Ns, theta_est_ABC);
load("Theta_true_values.mat") 
[Pv_true,t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_true);

t = t*1e9;

C = linspecer(3);

tiles = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(t,pow2db(mean(Pv_true,2)),'black')
hold on
plot(t,pow2db(mean(Pv_BSL,2)),'color',C(1,:),'LineWidth',2)
legend('True','BSL')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('[dB]')
box off
nexttile
plot(t,pow2db(mean(Pv_true,2)),'black')
hold on
plot(t,pow2db(mean(Pv_ABC_iter1,2)),'color',C(2,:),'LineWidth',2)
legend({sprintf('True'),sprintf('ABC\n1 iter')})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('[dB]')
box off
nexttile
plot(t,pow2db(mean(Pv_true,2)),'black')
hold on
plot(t,pow2db(mean(Pv_ABC_iter20,2)),'color',C(3,:),'LineWidth',2)
legend({sprintf('True'),sprintf('ABC\n20 iter')})
xtickformat('%d ns')
ylabel('[dB]')
box off
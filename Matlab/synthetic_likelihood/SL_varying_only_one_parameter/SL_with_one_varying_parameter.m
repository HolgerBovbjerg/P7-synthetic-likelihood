% clear all
%%
load('Prior_data_medium_prior_min_max_values.mat')
load('Theta_true_values.mat')
load('Stats_obs_true.mat')

N = 200; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
L = 100;

%%
steps = 700;
G0s = linspace(3.8595e-09,4.28e-9,steps);
Ts = linspace(7.5e-9,8e-9,steps);
sigmas = linspace(1e-5,2e-5,steps);
% lambdas = linspace(prior(3,1),prior(3,2),steps);
%%
loglikelihood = zeros(steps,1);

tic
for j = 1:steps
    %     lambda = lambdas(j);
    %         T = Ts(j);
    sigma = sigmas(j);
    %     G0 = G0s(j)
    parfor i = 1:L
        theta_prop = [theta_true(1) theta_true(2) theta_true(4) theta_true(4)];
%         theta_prop = [theta_true(1) theta_true(2) lambda theta_true(4)];
%         theta_prop = [T theta_true(2) theta_true(3) theta_true(4)];
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop);
        s_sim(i,:) = create_statistics(Pv, t);
    end
    loglikelihood(j) = (synth_loglikelihood(s_obs,s_sim));
    j
end
time = toc

%%
figure(1)
hold on
plot(G0s,exp(loglikelihood_G0_L10_N200),'g')
plot(G0s,exp(loglikelihood_G0_L50_N200),'m')
plot(G0s,exp(loglikelihood_G0_L100_N200),'k')
legend('G0\_L10\_N200','G0\_L50\_N200','G0\_L100\_N200')
ylim([0 16e5])
% xline(theta_true(2),'r')
%%
figure(2)
hold on
plot(G0s,exp(loglikelihood_G0_L50_N200),'g')
plot(G0s,exp(loglikelihood_G0_L50_N100),'m')
plot(G0s,exp(loglikelihood_G0_L50_N50),'k')
legend('G0\_L50\_N200','G0\_L50\_N100','G0\_L50\_N50')
% xline(theta_true(2),'r')
%%
figure(3)
hold on
plot(G0s,exp(loglikelihood_G0_L50_N200),'g')
plot(G0s,10*exp(loglikelihood_G0_L200_N50),'m')
plot(G0s,200*exp(loglikelihood_G0_L1000_N10),'k')
legend('200 turins og 50 summaries','50 turins og 200 summaries (times 10 for scaling reasons)')
%%
% figure(5)
% hold on
% plot(G0s,exp(loglikelihood_G0_L200_N10),'m')
% plot(G0s,exp(loglikelihood_G0_L100_N10),'k')
% ylim([0 800])
% legend('L200 N10','L10 N200')
%%
figure(5)
hold on
plot(Ts,exp(loglikelihood_T_L200_N10),'m')
plot(Ts,exp(loglikelihood_T_L10_N200),'k')
ylim([0 800])
legend('L200 N10','L10 N200')
%%
figure(6)
hold on
plot(sigmas,exp(loglikelihood_Sigma_L100_N10),'m')
plot(exp(loglikelihood_Sigma_L200_N10),'k')
legend('sigmaN\_L100\ N10','sigmaN\_L200\ N10')
%%
figure(7)
hold on
% plot(G0s,exp(loglikelihood_Sigma_L10_N200),'g')
% plot(G0s,exp(loglikelihood_Sigma_L50_N200),'m')
plot(G0s,exp(loglikelihood_Sigma_L100_N200),'k')
legend('Sigma\_L10\_N200','Sigma\_L50\_N200','Sigma\_L100\_N200')
% ylim([0 16e5])
% xline(theta_true(2),'r')
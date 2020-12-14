 clear all
%%
load('Prior_data_medium_prior_min_max_values.mat')
load('Theta_true_values.mat')


N = 300; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
L = 500;

[Pv, t] = sim_turin_matrix_gpu(1000, B, Ns, theta_true);
s_obs = create_statistics_new(Pv, t);

%%
steps = 35;
% G0s = linspace(3.8595e-09,4.28e-9,steps);
% Ts = linspace(7.5e-9,8e-9,steps);
lambdas = linspace(prior(3,1),prior(3,2),steps);
%%
loglikelihood = zeros(steps,1);
tic
for j = 1:steps
    lambda = lambdas(j);

    tic
    parfor i = 1:L
        theta_prop = [1.222097614611416e-08 2.677507662312809e-08 lambda 1.951110337740534e-05];
%         theta_prop = [theta_true(1) G0 theta_true(3) theta_true(4)];
                
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop);
        s_sim(i,:) = create_statistics_new(Pv, t);
    end
    loglikelihood(j) = (synth_loglikelihood(s_obs,s_sim));
    j
%     steptime = toc
    disp(['j = ',num2str(j),' step time = ', num2str(toc)])
end
% toc
% 
% %%
% figure(1)
% hold on
% plot(G0s,exp(loglikelihood_G0_L10_N200),'g')
% plot(G0s,exp(loglikelihood_G0_L50_N200),'m')
% plot(G0s,exp(loglikelihood_G0_L100_N200),'k')
% legend('G0_L10','G0_L50','G0_L100')
% ylim([0 16e5])
% % xline(theta_true(2),'r')
% %%
% figure(2)
% hold on
% plot(G0s,exp(loglikelihood_G0_L50_N200),'g')
% plot(G0s,5exp(loglikelihood_G0_L50_N100),'m')
% plot(G0s,10exp(loglikelihood_G0_L50_N50),'k')
% legend('G0_L50_N200','G0_L50_N100','G0_L50_N50')
% % xline(theta_true(2),'r')
% %%
% figure(3)
% hold on
% plot(G0s,exp(loglikelihood_G0_L50_N200),'g')
% plot(G0s,10exp(loglikelihood_G0_L200_N50),'m')
% plot(G0s,200exp(loglikelihood_G0_L1000_N10),'k')
% legend('200 turins og 50 summaries','50 turins og 200 summaries (times 10 for scaling reasons)')
% %%
% figure(4)
% hold on
% plot(Ts,exp(loglikelihood_T_L10_N200),'o')
% plot(Ts,exp(loglikelihood_T_L10_N200_sim2),'o')
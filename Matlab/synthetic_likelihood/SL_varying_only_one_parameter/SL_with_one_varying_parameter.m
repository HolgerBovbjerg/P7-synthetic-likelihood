% clear all
%%
load('Prior_data_medium_prior_min_max_values.mat')
load('Theta_true_values.mat')
load('S_obs_true.mat')

N = 10; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
L = 200;

%%
steps = 700;
G0s = linspace(3.8595e-09,4.28e-9,steps);
Ts = linspace(7.5e-9,8e-9,steps);
% lambdas = linspace(prior(3,1),prior(3,2),steps);
%%
loglikelihood = zeros(steps,1);


tic
for j = 1:steps
%     lambda = lambdas(j);
    T = Ts(j);
    parfor i = 1:L
%         theta_prop = [theta_true(1) theta_true(2) lambda theta_true(4)];
        theta_prop = [T theta_true(2) theta_true(3) theta_true(4)];
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop);
        s_sim(i,:) = create_statistics(Pv, t);
    end
    loglikelihood(j) = (synth_loglikelihood(s_obs,s_sim));
    j
end
toc

%%
figure(1)
hold on
plot(G0s,exp(loglikelihood_G0_L10_N200),'g')
plot(G0s,exp(loglikelihood_G0_L50_N200),'m')
plot(G0s,exp(loglikelihood_G0_L100_N200),'k')
legend('G0\_L10','G0\_L50','G0\_L100')
ylim([0 16e5])
% xline(theta_true(2),'r')
%%
figure(2)
hold on
plot(G0s,exp(loglikelihood_G0_L50_N200),'g')
plot(G0s,5*exp(loglikelihood_G0_L50_N100),'m')
plot(G0s,10*exp(loglikelihood_G0_L50_N50),'k')
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
figure(4)
hold on
plot(Ts,exp(loglikelihood_T_L10_N200),'o')
plot(Ts,exp(loglikelihood_T_L10_N200_sim2),'o')
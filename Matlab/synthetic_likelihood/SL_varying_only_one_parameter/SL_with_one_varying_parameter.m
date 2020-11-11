clear all
load('Prior_data_medium_prior_min_max_values.mat')
load('Theta_true_values.mat')
load('S_obs_true.mat')

N = 200; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
L = 10;

steps = 60;
G0s = linspace(prior(2,1),prior(2,2),steps); 
lambdas = linspace(prior(3,1),prior(3,2),steps);

loglikelihood = zeros(steps,1);

k = 0; 
for j = 1:steps
%     lambda = lambdas(j);
    G0 = G0s(j);
    parfor i = 1:L
%         theta_prop = [theta_true(1) theta_true(2) lambda theta_true(4)];
        theta_prop = [theta_true(1) G0 theta_true(3) theta_true(4)];
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop);
        s_sim(i,:) = create_statistics(Pv, t);
    end
    loglikelihood(j) = (synth_loglikelihood(s_obs,s_sim));
    k = k + 1
end

plot(G0s,exp(loglikelihood))

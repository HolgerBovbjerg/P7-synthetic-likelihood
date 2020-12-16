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
steps = 10;
varies = 10;
% Ts = linspace(prior(1,1),prior(1,2),varies);
G0s = linspace(prior(2,1),prior(2,2),varies);
% Ts = linspace(7.5e-9,8e-9,steps);
% lambdas = linspace(prior(3,1),prior(3,2),steps);
lambdas = linspace(prior(3,1),prior(3,2),steps);
%%
loglikelihood = zeros(varies,steps,1);
for g = 1:varies
    G0 = G0s(g);
    for j = 1:steps
        lambda = lambdas(j);
        
        tic
        parfor i = 1:L
            %             theta_prop = [1.222097614611416e-08 2.677507662312809e-08 lambda 1.951110337740534e-05];
            theta_prop = [T theta_true(2) lambda 1.951110337740534e-05];
            
            [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop);
            s_sim(i,:) = create_statistics_new(Pv, t);
        end
        loglikelihood(g,j) = (synth_loglikelihood(s_obs,s_sim));
        disp(['j = ',num2str(j),' step time = ', num2str(toc)])
    end
end

%%

ttt = tiledlayout(4,3);

for i = 1:10
    nexttile
    hold on
    plot(lambdas,loglikelihood(i,:))
end











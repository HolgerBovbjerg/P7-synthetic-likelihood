%%

L = 500;
s_obs = zeros(4,L);
mean_obs = zeros(4,L);
for i = 1:L
    i
    [Pv, t] = sim_turin_matrix_gpu(200, B, Ns, theta_true);
    s_obs(:,i) = create_statistics(Pv, t);
    mean_obs(:,i) = mean(s_obs(:,1:i),2);
end
%%
tiledlayout(4,1)
nexttile
plot(mean_obs(1,:))

nexttile
plot(mean_obs(2,:))

nexttile
plot(mean_obs(3,:))

nexttile
plot(mean_obs(4,:))

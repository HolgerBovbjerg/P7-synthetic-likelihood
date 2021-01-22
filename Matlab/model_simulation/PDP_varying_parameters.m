clear

%% clear
clear 

N = 300;% % Number of data sets to generate and average over
B = 4e9; % Bandwidth of signal: 4 GHz
steps = 5; % Number of different lambda values
Ns = 801; % Number of sample points in each data set
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
lambda = 1e9; % Arrival rate lambda 10e9 arrivals per second
sigma_N = sqrt(0.28e-9); % Noise standard deviation

%% Plot
tiledlayout(2,2)

C = linspecer(5);
% T
nexttile
Ts = [ 1e-9 4e-9 7.8e-9 10e-9 15e-9]; 
P_y_simulated = zeros(Ns,steps);
for i = 1:steps
    [P_y, t] = sim_turin_matrix_gpu(N,B,Ns,[Ts(i) G0 lambda sigma_N]); % N realisations of Turin model power delay profile
    P_y_simulated(:,i) = mean(P_y,2); 
end
t = t*1e9;
Tvals = ["1\cdot10^{-9}", "4\cdot10^{-9}", "7.8\cdot10^{-9}", "10\cdot10^{-9}", "15\cdot10^{-9}"];
for i = 1:steps
    plot(t,pow2db(P_y_simulated(:,i)),'DisplayName', "T= " + Tvals(i), 'color',C(i,:),'LineWidth',2)
    hold on
end
legend
xtickformat('%d ns')
ylabel('[dB]')
box off
% G0
nexttile
G0s = [db2pow(-74) db2pow(-80) db2pow(-84) db2pow(-90) db2pow(-94)]; 
P_y_simulated = zeros(Ns,steps);
for i = 1:steps
    [P_y, t] = sim_turin_matrix_gpu(N,B,Ns,[T G0s(i) lambda sigma_N]); % N realisations of Turin model power delay profile
    P_y_simulated(:,i) = mean(P_y,2); 
end
t = t*1e9;
G0vals = ["-74 [dB]", "-80 [dB]", "-84 [dB]", "-90 [dB]", "-94 [dB]"];
for i = 1:steps
    plot(t,pow2db(P_y_simulated(:,i)),'DisplayName', "G_0 = " + G0vals(i), 'color',C(i,:),'LineWidth',2)
    hold on
end
legend
xtickformat('%d ns')
ylabel('[dB]')
box off
% lambda
nexttile
lambdas = [5e6 1e7 1e8 1e9 4e9]; 
P_y_simulated = zeros(Ns,steps);
for i = 1:steps-1
    [P_y, t] = sim_turin_matrix_gpu(N,B,Ns,[T G0 lambdas(i) sigma_N]); % N realisations of Turin model power delay profile
    P_y_simulated(:,i) = mean(P_y,2); 
end
t = t*1e9;
lambdavals = ["5\cdot10^6", "10\cdot10^6", "100\cdot10^6", "1\cdot10^9", "4\cdot10^9"];
for i = 1:steps-1
    plot(t,pow2db(P_y_simulated(:,i)),'DisplayName', "\lambda = " + lambdavals(i), 'color',C(i,:),'LineWidth',2)
    hold on
end
legend
xtickformat('%d ns')
ylabel('[dB]')
box off
% sigma_N
nexttile

sigma_Ns = [sqrt(0.28e-7) sqrt(0.28e-8) sqrt(0.28e-9) sqrt(0.28e-10) sqrt(0.28e-11)];
P_y_simulated = zeros(Ns,steps);
for i = 1:steps
    [P_y, t] = sim_turin_matrix_gpu(N,B,Ns,[T G0 lambda sigma_Ns(i)]); % N realisations of Turin model power delay profile
    P_y_simulated(:,i) = mean(P_y,2); 
end
t = t*1e9;
sigma_N_square_vals = ["0.28\cdot10^{-7}", "0.28\cdot10^{-8}", "0.28\cdot10^{-9}", "0.28\cdot10^{-10}", "0.28\cdot10^{-11}"];
for i = 1:steps
    plot(t,pow2db(P_y_simulated(:,i)),'DisplayName', "\sigma_N^2 = " + sigma_N_square_vals(i), 'color',C(i,:),'LineWidth',2)
    hold on
end
legend
xtickformat('%d ns')
ylabel('[dB]')
box off
%% Distribution of summary statistics
clear
T = 7.8e-9; % Reverberation time: 7.8 ns
G0 = db2pow(-83.9); % Reverberation gain converted from dB to power
sigma_N = sqrt(28e-9); % Noise variance
B = 4e9; % Bandwidth of signal: 4 GHz
Ns = 801; % Number of frequency samples in transfer function

% Time delay tau is a possion arrival process with mean delay lambda
lambda = 10e9; % randomly chosen arrival rate lambda 10e9 arrivals per second

Nr = 100;    % Number of Turin realizations pr summary statistic vector
Nl = 1000;   % Number of different summary statistic vectors generated

ss1 = zeros(3,Nl,8); % Matrix for summary statistics

tic
temp = create_statistics(Nr, T, G0, lambda, sigma_N, B, Ns);
t_per_stat = toc;


est_time = Nl*t_per_stat

%%

T_min = 7.8e-11; T_max = 7.8e-7;
T = T_min + (T_max-T_min).*rand(Nl,1);

G0_min = db2pow(-110); G0_max = db2pow(-60);
G0 = G0_min + (G0_max-G0_min).*rand(Nl,1);
G0 = sort(G0);
sigma_N_min = sqrt(5e-10); sigma_N_max = sqrt(5e-9);
sigma_N = sigma_N_min + (sigma_N_max-sigma_N_min).*rand(Nl,1);

lambda_min = 50e6; lambda_max = 50e7;
lambda = lambda_min + (lambda_max-lambda_min).*rand(Nl,1);

%%
tic
for j = 1:3
    for i = 1:Nl    
        ss1(j,i,:) = create_statistics(Nr, T(j), G0(j), lambda(j), sigma_N(j), B, Ns);
    end
end 
toc
%%
bins = 35;

for j = 1:3
    figure(j)
    plotstatistics =  tiledlayout(2,4);
    title(plotstatistics,"Histograms of summary statistics, from "+ Nl +" summary statistic vectors, with random input parameters") 

    for k = 1:4
        nexttile
        hold on
        pd = fitdist(ss1(j,:,k)','Normal');
        x_vals = linspace(min(ss1(j,:,k)),max(ss1(j,:,k)),bins);
        mu = mean(pd);
        sigma = sqrt(var(pd));
        histogram(ss1(j,:,k),'Normalization','pdf','Binedges',x_vals,'FaceColor','b')
        p = plot(x_vals,pdf(pd,x_vals),'r');
        legend(p,"\mu = "+num2str(mu,2)+", \sigma = "+num2str(sigma,2))
        title("Histogram of mean of temporal moment "+(k-1))
    end
    
    for k = 5:8
        nexttile
        hold on
        pd = fitdist(ss1(j,:,k)','Normal');
        x_vals = linspace(min(ss1(j,:,k)),max(ss1(j,:,k)),bins);
        mu = mean(pd);
        sigma = sqrt(var(pd));
        histogram(ss1(j,:,k),'Normalization','pdf','Binedges',x_vals,'FaceColor','b')
        p = plot(x_vals,pdf(pd,x_vals),'r');
        legend(p,"\mu = "+num2str(mu,2)+", \sigma = "+num2str(sigma,2))
        title("Histogram of variance of temporal moment "+(k-5))
    end   
    
%     nexttile
%     hold on
%     pd = fitdist(ss1(j,:,2)','Normal');
%     plot(linspace(min(ss1(j,:,2)),max(ss1(j,:,2)),bins),pdf(pd,linspace(min(ss1(j,:,2)),max(ss1(j,:,2)),bins)))
%     histogram(ss1(j,:,2),'Normalization','pdf','Binedges',linspace(min(ss1(j,:,2)),max(ss1(j,:,2)),bins))
%     title("Histogram of mean of 1st temporal moment")
% 
%     nexttile
%     hold on
%     pd = fitdist(ss1(j,:,3)','Normal');
%     plot(linspace(min(ss1(j,:,3)),max(ss1(j,:,3)),bins),pdf(pd,linspace(min(ss1(j,:,3)),max(ss1(j,:,3)),bins)))
%     histogram(ss1(j,:,3),'Normalization','pdf','Binedges',linspace(min(ss1(j,:,3)),max(ss1(j,:,3)),bins))
%     title("Histogram of mean of 2nd temporal moment")
% 
%     nexttile
%     hold on
%     pd = fitdist(ss1(j,:,4)','Normal');
%     plot(linspace(min(ss1(j,:,4)),max(ss1(j,:,4)),bins),pdf(pd,linspace(min(ss1(j,:,4)),max(ss1(j,:,4)),bins)))
%     histogram(ss1(j,:,4),'Normalization','pdf','Binedges',linspace(min(ss1(j,:,4)),max(ss1(j,:,4)),bins))
%     title("Histogram of mean of 3rd temporal moment")
% 
%     nexttile
%     hold on
%     pd = fitdist(ss1(j,:,5)','Normal');
%     plot(linspace(min(ss1(j,:,5)),max(ss1(j,:,5)),bins),pdf(pd,linspace(min(ss1(j,:,5)),max(ss1(j,:,5)),bins)))
%     histogram(ss1(j,:,5),'Normalization','pdf','Binedges',linspace(min(ss1(j,:,5)),max(ss1(j,:,5)),bins))
%     title("Histogram of variance of 0th temporal moment")
% 
%     nexttile
%     hold on
%     pd = fitdist(ss1(j,:,6)','Normal');
%     plot(linspace(min(ss1(j,:,6)),max(ss1(j,:,6)),bins),pdf(pd,linspace(min(ss1(j,:,6)),max(ss1(j,:,6)),bins)))
%     histogram(ss1(j,:,6),'Normalization','pdf','Binedges',linspace(min(ss1(j,:,6)),max(ss1(j,:,6)),bins))
%     title("Histogram of variance of 1st temporal moment")
% 
%     nexttile
%     hold on
%     pd = fitdist(ss1(j,:,7)','Normal');
%     plot(linspace(min(ss1(j,:,7)),max(ss1(j,:,7)),bins),pdf(pd,linspace(min(ss1(j,:,7)),max(ss1(j,:,7)),bins)))
%     histogram(ss1(j,:,7),'Normalization','pdf','Binedges',linspace(min(ss1(j,:,7)),max(ss1(j,:,7)),bins))
%     title("Histogram of variance of 2nd temporal moment")
% 
%     nexttile
%     hold on
%     pd = fitdist(ss1(j,:,8)','Normal');
%     plot(linspace(min(ss1(j,:,8)),max(ss1(j,:,8)),bins),pdf(pd,linspace(min(ss1(j,:,8)),max(ss1(j,:,8)),bins)))
%     histogram(ss1(j,:,8),'Normalization','pdf','Binedges',linspace(min(ss1(j,:,8)),max(ss1(j,:,8)),bins))
%     title("Histogram of variance of 3rd temporal moment")
end
%%
pd = fitdist(ss1(j,:,1)','Normal');
figure(1)
hold on
% histogram(ss1(j,:,1),'Binedges',linspace(min(ss1(j,:,1)),max(ss1(j,:,1)),bins))
histogram(ss1(j,:,1),'Normalization','pdf','Binedges',linspace(min(ss1(j,:,1)),max(ss1(j,:,1)),bins))
plot(linspace(min(ss1(j,:,1)),max(ss1(j,:,1)),bins),pdf(pd,linspace(min(ss1(j,:,1)),max(ss1(j,:,1)),bins)))
title("Histogram of mean of 0th temporal moment")
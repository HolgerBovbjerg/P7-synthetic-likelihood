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
Nl = 600;   % Number of different summary statistic vectors generated

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
for j = 1:6
    for i = 1:Nl    
        ss1(j,i,:) = create_statistics(Nr, T(j), G0(j), lambda(j), sigma_N(j), B, Ns);
    end
end 
toc
%%
bins = 25;

for j = 1:6
    figure(j)
    plotstatistics(j) =  tiledlayout(2,4);
    title(plotstatistics,"Histograms of summary statistics, from "+ Nl +" summary statistic vectors, with random input parameters") 

    for k = 1:4
        nexttile
        hold on
        pd_normal = fitdist(ss1(j,:,k)','Normal');
        pd_lognormal = fitdist(ss1(j,:,k)','Lognormal');
        x_vals = linspace(min(ss1(j,:,k)),max(ss1(j,:,k)),bins);
        histogram(ss1(j,:,k),'Normalization','pdf','Binedges',x_vals,'FaceColor','b')
        p_normal = plot(x_vals,pdf(pd_normal,x_vals),'r', 'Linewidth',3);
        p_lognormal = plot(x_vals,pdf(pd_lognormal,x_vals),'g', 'Linewidth',3);
        legend([p_normal p_lognormal],{'Normal','Lognormal'})
%         legend([p_normal p_lognormal],{"\mu = "+num2str(pd_normal.mu,2)+", \sigma = "+num2str(pd_normal.sigma,2),'Lognormal'})
        title("     Histogram of mean of temporal moment "+(k-1))
    end
    
    for k = 5:8
        nexttile
        hold on
        pd_normal = fitdist(ss1(j,:,k)','Normal');
        pd_lognormal = fitdist(ss1(j,:,k)','Lognormal');
        x_vals = linspace(min(ss1(j,:,k)),max(ss1(j,:,k)),bins);
        histogram(ss1(j,:,k),'Normalization','pdf','Binedges',x_vals,'FaceColor','b')
        p_normal = plot(x_vals,pdf(pd_normal,x_vals),'r', 'Linewidth',3);
        p_lognormal = plot(x_vals,pdf(pd_lognormal,x_vals),'g', 'Linewidth',3);
        legend([p_normal p_lognormal],{'Normal','Lognormal'})
%         legend([p_normal p_lognormal],{"\mu = "+num2str(pd_normal.mu,2)+", \sigma = "+num2str(pd_normal.sigma,2),'Lognormal'})
        title("          Histogram of variance of temporal moment "+(k-5))
    end   
end
%%
% for j = 1:6
%     exportgraphics(plotstatistics(j),"Distribution_of_summary_statistics_with_fixed_parameter_set_"+num2str(j)+".pdf",'ContentType','vector')
% end
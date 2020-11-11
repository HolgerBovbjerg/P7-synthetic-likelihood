%%

figure(5)
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
blocks = 1;
thetas_plot = thetas(:,1:j);
tt = tiledlayout(4,blocks)
title(tt,"Density plots ")
len = length(thetas_plot(1,1:j));

linesize = 1;
kslinesizefactor = 1.5;
fontsize =16;

T_true = theta_true(1);
T_start = thetas(1,1);

G0_true = theta_true(2);
G0_start = thetas(2,1);

lambda_true = theta_true(3);
lambda_start = thetas(3,1);

sigma_N_true = theta_true(4);
sigma_start = thetas(4,1);

%prior = prior_very_small;

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(1,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(1,1:len/blocks*(kk)),linspace(prior(1,1), prior(1,2),1000));
    plot(x,f,'LineWidth',linesize,'Color','black')
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    %plot(MAPx,MAPy,'Color','red')
    xline(x(index),'Color','red')
    xline(T_true,'--','LineWidth',linesize)
    %xline(T_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(1,1) prior(1,2)])
    title("$T$",'Fontsize',fontsize)
end

legend( "Approx. posterior", "MAP estimate","True value")

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(2,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(2,1:len/blocks*(kk)),linspace(prior(2,1), prior(2,2),1000));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    xline(x(index),'Color','red')
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(G0_true,'--','LineWidth',linesize)
    %xline(G0_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(2,1) prior(2,2)])
    title("$G_0$",'Fontsize',fontsize)
end

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(3,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(3,1:len/blocks*(kk)),linspace(prior(3,1), prior(3,2),1000));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    xline(x(index),'Color','red')
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(lambda_true,'--','LineWidth',linesize)
    %xline(lambda_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(3,1) prior(3,2)])
    title("$\lambda$",'Fontsize',fontsize)
end

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(4,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(4,1:len/blocks*(kk)),linspace(prior(4,1), prior(4,2),1000));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    xline(x(index),'Color','red')
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(sigma_N_true,'--','LineWidth',linesize)
    %xline(sigma_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(4,1) prior(4,2)])
    xlabel([num2str(round((len/blocks)*kk)),' steps'])
    title("$\sigma_N$",'Fontsize',fontsize)
end

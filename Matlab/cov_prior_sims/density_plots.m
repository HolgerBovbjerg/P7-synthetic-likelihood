%%

figure(3)
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
blocks = 1;

thetas_plot = thetas(:,1:j);

tt = tiledlayout(2,2)
% title(tt,"Density plots - small covariance")
    
len = length(thetas_plot(1,1:j));

linesize = 1;
kslinesizefactor = 1.5;
fontsize =16;

burnin = 200;

T_true = theta_true(1);
T_start = thetas(1,1);

G0_true = theta_true(2);
G0_start = thetas(2,1);

lambda_true = theta_true(3);
lambda_start = thetas(3,1);

sigma_N_true = theta_true(4);
sigma_start = thetas(4,1);

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(1,1+len/blocks*(kk-1):len/blocks*(kk)));

    [f x bw] = ksdensity(thetas_plot(1,burnin:len/blocks*(kk)),linspace(prior(1,1), prior(1,2),100000));
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(mean(thetas(1,burnin:j)))
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    %plot(MAPx,MAPy,'Color','red')
    xline(x(index),'Color','red')
    xline(T_true,'--','LineWidth',linesize)
    %xline(T_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(1,1) prior(1,2)])
    xticks(prior(1,1):prior(1,2)/5:prior(1,2))
    title("$T$",'Fontsize',fontsize)
    set(gca,'Yticklabel',[])
%     set(tt,'ytick',[])

end

legend( "Approx. posterior",'Mean', "MAP estimate","True value")

thetas_plot(2,:) = pow2db(thetas_plot(2,:));

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(2,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(2,burnin:len/blocks*(kk)),pow2db(linspace(prior(2,1), prior(2,2),100000)));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    xline(x(index),'Color','red')
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(mean(thetas_plot(2,burnin:j-1)))
    xline(pow2db(G0_true),'--','LineWidth',linesize)
    %xline(G0_start,'LineWidth',linesize,'Color','blue')
    xlim([pow2db(prior(2,1)) pow2db(prior(2,2))])
    xticks(pow2db(prior(2,1)):2:pow2db(prior(2,2)))
    set(gca,'Yticklabel',[])
    title("$G_0$",'Fontsize',fontsize)
end

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(3,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(3,burnin:len/blocks*(kk)),linspace(prior(3,1), prior(3,2),100000));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    xline(x(index),'Color','red')
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(mean(thetas(3,burnin:j)))
    xline(lambda_true,'--','LineWidth',linesize)
    %xline(lambda_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(3,1) prior(3,2)])
    xticks(prior(3,1):prior(3,2)/5:prior(3,2))
    set(gca,'Yticklabel',[])
    title("$\lambda$",'Fontsize',fontsize)
end

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(4,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(4,burnin:len/blocks*(kk)),linspace(prior(4,1), prior(4,2),100000));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    xline(x(index),'Color','red')
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(mean(thetas(4,burnin:j)))
    xline(sigma_N_true,'--','LineWidth',linesize)
    %xline(sigma_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(4,1) prior(4,2)])
    xticks(prior(4,1):prior(4,2)/5:prior(4,2))
    set(gca,'Yticklabel',[])
    title("$\sigma_N$",'Fontsize',fontsize)
end

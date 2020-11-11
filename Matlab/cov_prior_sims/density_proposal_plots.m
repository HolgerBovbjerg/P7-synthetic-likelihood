%%
scale = 1;
normalT = makedist("Normal","mu",theta_true(1),"sigma",sqrt(covariance(1,1)/scale))
normalG0 = makedist("Normal","mu",theta_true(2),"sigma",sqrt(covariance(2,2)/scale))
normallambda = makedist("Normal","mu",theta_true(3),"sigma",sqrt(covariance(3,3)/scale))
normalsigma = makedist("Normal","mu",theta_true(4),"sigma",sqrt(covariance(4,4)/scale))


yT = linspace(prior(1,1),prior(1,2),500);
yG0 = linspace(prior(2,1),prior(2,2),500);
ylambda = linspace(prior(3,1),prior(3,2),500);
ysigma = linspace(prior(4,1),prior(4,2),500);

figure(5)
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
blocks = 1;
thetas_plot = thetas(:,1:j);
tt = tiledlayout(4,blocks)
title(tt,"Density plots - small prior")
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
    plot(yT,pdf(normalT,yT),'b--')
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

legend( "Approx. posterior", "Proposal distribution", "MAP estimate","True value")

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
    plot(yG0,pdf(normalG0,yG0),'b--')
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
    plot(ylambda,pdf(normallambda,ylambda),'b--')
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
    plot(ysigma,pdf(normalsigma,ysigma),'b--')
    xline(sigma_N_true,'--','LineWidth',linesize)
    %xline(sigma_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(4,1) prior(4,2)])
    xlabel([num2str(round((len/blocks)*kk)),' steps'])
    title("$\sigma_N$",'Fontsize',fontsize)
end


%%

ttt = tiledlayout(4,1);
title(ttt,"Accepted parameters in Markov Chain using medium prior",'Fontsize',24)
subtitle(ttt,['Number of steps: ',num2str(j),'  Acceptance rate: ',num2str(100*accept/j),' %'])
nexttile
plot(thetas(1,1:j),'o')
ylim([prior(1,1) prior(1,2)])
yline(theta_true(1))
ylabel("Magnitude")
title("T")

nexttile
plot(pow2db(thetas(2,1:j)),'o')
ylim([pow2db(prior(2,1)) pow2db(prior(2,2))])
yline(pow2db(theta_true(2)))
ylabel("dB")
title("G0")

nexttile
plot(thetas(3,1:j),'o')
ylim([prior(3,1) prior(3,2)])
yline(theta_true(3))
ylabel("Magnitude")
title("$\lambda$")

nexttile
plot(thetas(4,1:j),'o')
ylim([prior(4,1) prior(4,2)])
yline(theta_true(4))
ylabel("Magnitude")
title("$\sigma$")

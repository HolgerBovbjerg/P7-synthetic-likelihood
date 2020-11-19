%%

% figure(3)
set(0,'defaultTextInterpreter','tex');
set(groot, 'defaultLegendInterpreter','tex');
blocks = 1;

thetas_plot = thetas(:,1:j);

tt = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
% subtitle(tt,"Density plots - small covariance")
% tt.Units = 'inches';
% tt.OuterPosition = [0 0 3.5 3.5];    
len = length(thetas_plot(1,1:j));

linesize = 2;
fontsize = 8;

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
    %[f, x] = ksdensity(thetas(1,1+len/blocks*(kk-1):len/blocks*(kk)));

    [f, x] = ksdensity(thetas_plot(1,burnin:len/blocks*(kk)),linspace(prior(1,1), prior(1,2),100000));
    area(x,f,'LineWidth',1,'FaceColor','#bbbbbb')
    xline(mean(thetas(1,burnin:j)),'r','LineWidth',linesize)
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    xline(T_true,'--g','LineWidth',linesize)
    xlim([prior(1,1) prior(1,2)])
    xticks([prior(1,1) 4.5e-9 8e-9 11.5e-9 prior(1,2)])
    subtitle('T','Fontsize',fontsize)
    set(gca,'Yticklabel',[])
    set(gca,'ytick',[])
end

legend( "Approx. posterior",'MMSEE',"True value")

thetas_plot(2,:) = pow2db(thetas_plot(2,:));

for kk = 1:blocks
    nexttile
    hold on
    %[f, x] = ksdensity(thetas(2,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f, x] = ksdensity(thetas_plot(2,burnin:len/blocks*(kk)),pow2db(linspace(prior(2,1), prior(2,2),100000)));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
%     xline(x(index),'Color','red')
    area(x,f,'LineWidth',1,'FaceColor','#bbbbbb')
    xline(mean(thetas_plot(2,burnin:j-1)),'r','LineWidth',linesize)
    xline(pow2db(G0_true),'--g','LineWidth',linesize)
    %xline(G0_start,'LineWidth',linesize,'Color','blue')
    xlim([pow2db(prior(2,1)) pow2db(prior(2,2))])
    xticks([pow2db(prior(2,1)) -89 -84 -79 pow2db(prior(2,2))])
    set(gca,'Yticklabel',[])
    set(gca,'ytick',[])
    subtitle('G_0','Fontsize',fontsize)
end

for kk = 1:blocks
    nexttile
    hold on
    %[f, x] = ksdensity(thetas(3,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f, x] = ksdensity(thetas_plot(3,burnin:len/blocks*(kk)),linspace(prior(3,1), prior(3,2),100000));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
%     xline(x(index),'Color','red')
    area(x,f,'LineWidth',1,'FaceColor','#bbbbbb')
    xline(mean(thetas(3,burnin:j)),'r','LineWidth',linesize)
    xline(lambda_true,'--g','LineWidth',linesize)
    %xline(lambda_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(3,1) prior(3,2)])
    xticks([prior(3,1) 1e9 2e9 3e9 prior(3,2)])
    set(gca,'Yticklabel',[])
    set(gca,'ytick',[])
    subtitle('\lambda','Fontsize',fontsize)
end

for kk = 1:blocks
    nexttile
    hold on
    %[f, x] = ksdensity(thetas(4,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f, x] = ksdensity(thetas_plot(4,burnin:len/blocks*(kk)),linspace(prior(4,1), prior(4,2),100000));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
%     xline(x(index)^2,'Color','red')
    area(x.^2,f,'LineWidth',1,'FaceColor','#bbbbbb')
    xline(mean(thetas(4,burnin:j))^2,'r','LineWidth',linesize)
    xline(sigma_N_true^2,'--g','LineWidth',linesize)
    %xline(sigma_start,'LineWidth',linesize,'Color','blue')
    xlim([prior(4,1)^2 prior(4,2)^2])
    xticks([prior(4,1)^2 0.7e-9 1.4e-9 2.1e-9 prior(4,2)^2])
    set(gca,'Yticklabel',[])
    set(gca,'ytick',[])
    subtitle('\sigma_N^2','Fontsize',fontsize)
end

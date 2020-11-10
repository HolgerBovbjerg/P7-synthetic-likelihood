%%

figure(5)
blocks = 5;
thetas_plot = thetas(:,1:j);
tt = tiledlayout(4,blocks)
title(tt,"Density plots ")
len = length(thetas_plot(1,1:j));

linesize = 1;
kslinesizefactor = 1.5;


% T_true = 7.8e-9;
T_start = thetas(1,1);

% G0_true = 4.0738e-09;
G0_start = thetas(2,1);

% lambda_true = 10e8;
lambda_start = thetas(3,1);

% sigma_true = 1.6733e-5;
sigma_start = thetas(4,1);

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(1,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(1,1:len/blocks*(kk)));
    plot(x,f,'LineWidth',linesize,'Color','black')
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    plot(MAPx,MAPy,'Color','red')
    xline(T_true,'--','LineWidth',linesize)
    xline(T_start,'LineWidth',linesize,'Color','blue')
    xline(mean(thetas(1,:)))
   % xlim([2e-9 9e-9])
    title("T")
end

legend( "Approx. posterior", "MAP estimate","True value", "Start value")

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(2,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(2,1:len/blocks*(kk)));
    [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    plot(MAPx,MAPy,'Color','red')
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(G0_true,'--','LineWidth',linesize)
    xline(G0_start,'LineWidth',linesize,'Color','blue')
    xline(mean(thetas(2,:)))
   % xlim([db2pow(-85) db2pow(-82)])
    title("G_0")
end

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(3,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(3,1:len/blocks*(kk)));
        [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    plot(MAPx,MAPy,'Color','red')
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(lambda_true,'--','LineWidth',linesize)
    xline(lambda_start,'LineWidth',linesize,'Color','blue')
    xline(mean(thetas(3,:)))
   % xlim([5e8 15e8])
    title("\lambda")
end

for kk = 1:blocks
    nexttile
    hold on
    %[f x] = ksdensity(thetas(4,1+len/blocks*(kk-1):len/blocks*(kk)));
    [f x] = ksdensity(thetas_plot(4,1:len/blocks*(kk)));
        [MAPest, index] = max(f);
    MAPx = ones(100,1)*x(index);
    MAPy = linspace(0,MAPest,100);
    plot(MAPx,MAPy,'Color','red')
    plot(x,f,'LineWidth',linesize,'Color','black')
    xline(sigma_N_true,'--','LineWidth',linesize)
    xline(sigma_start,'LineWidth',linesize,'Color','blue')
    xline(mean(thetas(4,:)))
  %  xlim([sqrt(0.1e-9) sqrt(0.5e-9)])
    xlabel([num2str(round((len/blocks)*kk)),' steps'])
    title("\sigma_N")
end
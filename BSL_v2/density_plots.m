%%

figure(1)
blocks = 5;
tiledlayout(4,blocks)
len = length(thetas(1,1:end));

linesize = 2;
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
    [f x] = ksdensity(thetas(1,1+len/blocks*(kk-1):len/blocks*(kk)));
    plot(x,f,'LineWidth',linesize*kslinesizefactor)
    xline(T_true,'--','LineWidth',linesize)
    xline(T_start,'LineWidth',linesize,'Color','blue')
    %xlim([4e-9 10e-9])
    title("T")
end

for kk = 1:blocks
    nexttile
    [f x] = ksdensity(thetas(2,1+len/blocks*(kk-1):len/blocks*(kk)));
    plot(x,f,'LineWidth',linesize*kslinesizefactor)
    xline(G0_true,'--','LineWidth',linesize)
    xline(G0_start,'LineWidth',linesize,'Color','blue')
   % xlim([2e-9 8e-9])
    title("G0")
end

for kk = 1:blocks
    nexttile
    [f x] = ksdensity(thetas(3,1+len/blocks*(kk-1):len/blocks*(kk)));
    plot(x,f,'LineWidth',linesize*kslinesizefactor)
    xline(lambda_true,'--','LineWidth',linesize)
    xline(lambda_start,'LineWidth',linesize,'Color','blue')
   % xlim([3e8 16e8])
    title("\lambda")
end

for kk = 1:blocks
    nexttile
    [f x] = ksdensity(thetas(1,1+len/blocks*(kk-1):len/blocks*(kk)));
    plot(x,f,'LineWidth',linesize*kslinesizefactor)
    xline(sigma_true,'--','LineWidth',linesize)
    xline(sigma_start,'LineWidth',linesize,'Color','blue')
   % xlim([0.2e-5 2e-5])
    title("\sigma")
end

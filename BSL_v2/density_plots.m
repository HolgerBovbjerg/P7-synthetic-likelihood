%%

figure(1)
blocks = 5;
tiledlayout(4,blocks)
len = length(thetas(1,1:1000));

linesize = 2;
kslinesizefactor = 1.5;


T_true = 7.8e-9;
T_start = T_true/1.2;

G0_true = 4.0738e-09;
G0_start = G0_true/1.2;

lambda_true = 10e8;
lambda_start = lambda_true*1.2;

sigma_true = 1.6733e-5;
sigma_start = sigma_true*1.2;

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

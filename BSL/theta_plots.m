  %%
      figure(1)
    tiledlayout(4,1)
    nexttile
    
    plot(1:length(values),values(1,:),"o")
    yline(7.8e-9)
   % ylim([6e-9 9e-9])
    title("T")
    
    nexttile
    
    plot(1:length(values),values(2,:),"o")
    yline(db2pow(-83.9))
   % ylim([3e-9 6e-9])
    title("G0")
    
    nexttile
    
    plot(1:length(values),values(3,:),"o")
    yline(10e8)
   % ylim([0.5e+9 1.5e+9])
    title("\lambda")
    
    nexttile
    
    plot(1:length(values),values(4,:),"o")
    yline(sqrt(0.28e-9))
    %ylim([1.0e-5 2.5e-5])
    title("\sigma")

    
    %%   Plotting independent proposal distribution
    
load("thetas_independent.mat")
indvalues = values;
clear values;
      figure(4)
      linesize = 2;
      kslinesizefactor = 1.5;
    tiledlayout(4,8)
    
Tstart = 5.055e-9;
Ttrue = 7.8e-9;
G0start = 6.056e-09; % Reverberation gain converted from dB to power
G0true = 4.0738e-09;

nexttile
[f x] = ksdensity(indvalues(1,2:1000));
plot(x,f,'LineWidth',linesize*kslinesizefactor)
xline(Ttrue,'--','LineWidth',linesize)
xline(Tstart,'LineWidth',linesize,'Color','blue')
xlim([4e-9 10e-9])
title("T")

for kk = 1:7
    nexttile
    [f x] = ksdensity(indvalues(1,1004+(kk-1)*1500:1004+kk*1500));
    plot(x,f,'LineWidth',linesize*kslinesizefactor)
    xline(Ttrue,'--','LineWidth',linesize)
    xline(Tstart,'LineWidth',linesize,'Color','blue')
    xlim([4e-9 10e-9])
    title("G0")
end






%--------------------------------------------------------
G0start = 6.056e-09; % Reverberation gain converted from dB to power
G0true = 4.0738e-09;

nexttile
[f x] = ksdensity(indvalues(2,2:1000));
plot(x,f,'LineWidth',linesize*kslinesizefactor)
xline(G0true,'--','LineWidth',linesize)
xline(G0start,'LineWidth',linesize,'Color','blue')
xlim([2e-9 8e-9])
title("\lambda")

for kk = 1:7
    nexttile
    [f x] = ksdensity(indvalues(2,1004+(kk-1)*1500:1004+kk*1500));
    plot(x,f,'LineWidth',linesize*kslinesizefactor)
    xline(G0true,'--','LineWidth',linesize)
    xline(G0start,'LineWidth',linesize,'Color','blue')
    xlim([2e-9 8e-9])
    title("G0")
end

%--------------------------------------------------------
lambdatrue = 10e8;
lambdastart = 0.599878e+09;

nexttile
[f x] = ksdensity(indvalues(3,2:1000));
plot(x,f,'LineWidth',linesize*kslinesizefactor)
xline(lambdatrue,'--','LineWidth',linesize)
xline(lambdastart,'LineWidth',linesize,'Color','blue')
xlim([3e8 16e8])
title("\lambda")

for kk = 1:7
    nexttile
    [f x] = ksdensity(indvalues(3,1004+(kk-1)*1500:1004+kk*1500));
    plot(x,f,'LineWidth',linesize*kslinesizefactor)
    xline(lambdatrue,'--','LineWidth',linesize)
    xline(lambdastart,'LineWidth',linesize,'Color','blue')
    xlim([3e8 16e8])
    title("\lambda")
end


%------------------------------------------------
sigmatrue = 1.6733e-05;
sigmastart = 1.073155e-05;

nexttile
[f x] = ksdensity(indvalues(4,2:1000));
plot(x,f,'LineWidth',linesize*kslinesizefactor)
xline(sigmatrue,'--','LineWidth',linesize)
xline(sigmastart,'LineWidth',linesize,'Color','blue')
xlim([0.2e-5 2e-5])
title("\sigma")

for kk = 1:7
    nexttile
    [f x] = ksdensity(indvalues(4,1004+(kk-1)*1500:1004+kk*1500));
    plot(x,f,'LineWidth',linesize*kslinesizefactor)
    xline(sigmatrue,'--','LineWidth',linesize)
    xline(sigmastart,'LineWidth',linesize,'Color','blue')
    xlim([0.2e-5 2e-5])
    title("\sigma")
end


%%   Plotting dependent proposal distribution
load("thetas_mvnrnd_big.mat")
mvnvalues = values;
clear values;
      figure(5)
      linesize = 2;
      kslinesizefactor = 1.5;
    tiledlayout(4,8)
    
Tstart = 5.055e-9;
Ttrue = 7.8e-9;
G0start = 6.056e-09; % Reverberation gain converted from dB to power
G0true = 4.0738e-09;

nexttile
[f x] = ksdensity(mvnvalues(1,2:1500));
plot(x,f,'LineWidth',linesize*kslinesizefactor,'Color','red')
xline(Ttrue,'--','LineWidth',linesize)
xline(Tstart,'LineWidth',linesize,'Color','blue')
xlim([4e-9 10e-9])
title("T")

for kk = 1:7
    nexttile
    [f x] = ksdensity(mvnvalues(1,1502+(kk-1)*1500:1503+kk*1497));
    plot(x,f,'LineWidth',linesize*kslinesizefactor,'Color','red')
    xline(Ttrue,'--','LineWidth',linesize)
    xline(Tstart,'LineWidth',linesize,'Color','blue')
    xlim([4e-9 10e-9])
    title("G0")
end






%--------------------------------------------------------
G0start = 6.056e-09; % Reverberation gain converted from dB to power
G0true = 4.0738e-09;

nexttile
[f x] = ksdensity(mvnvalues(2,2:1500));
plot(x,f,'LineWidth',linesize*kslinesizefactor,'Color','red')
xline(G0true,'--','LineWidth',linesize)
xline(G0start,'LineWidth',linesize,'Color','blue')
xlim([2e-9 8e-9])
title("\lambda")

for kk = 1:7
    nexttile
    [f x] = ksdensity(mvnvalues(2,1502+(kk-1)*1500:1503+kk*1497));
    plot(x,f,'LineWidth',linesize*kslinesizefactor,'Color','red')
    xline(G0true,'--','LineWidth',linesize)
    xline(G0start,'LineWidth',linesize,'Color','blue')
   % xlim([2e-9 8e-9])
    title("G0")
end

%--------------------------------------------------------
lambdatrue = 10e8;
lambdastart = 0.599878e+09;

nexttile
[f x] = ksdensity(mvnvalues(3,2:1500));
plot(x,f,'LineWidth',linesize*kslinesizefactor,'Color','red')
xline(lambdatrue,'--','LineWidth',linesize)
xline(lambdastart,'LineWidth',linesize,'Color','blue')
%xlim([3e8 16e8])
title("\lambda")

for kk = 1:7
    nexttile
    [f x] = ksdensity(mvnvalues(3,2:1500));
    plot(x,f,'LineWidth',linesize*kslinesizefactor,'Color','red')
    xline(lambdatrue,'--','LineWidth',linesize)
    xline(lambdastart,'LineWidth',linesize,'Color','blue')
  %  xlim([3e8 16e8])
    title("\lambda")
end


%------------------------------------------------
sigmatrue = 1.6733e-05;
sigmastart = 1.073155e-05;

nexttile
[f x] = ksdensity(abs(mvnvalues(4,2:1500)));
plot(x,f,'LineWidth',linesize*kslinesizefactor,'Color','red')
xline(sigmatrue,'--','LineWidth',linesize)
xline(sigmastart,'LineWidth',linesize,'Color','blue')
%xlim([0.2e-5 2e-5])
title("\sigma")

for kk = 1:7
    nexttile
    [f x] = ksdensity(mvnvalues(4,1502+(kk-1)*1500:1503+kk*1497));  
    plot(x,f,'LineWidth',linesize*kslinesizefactor,'Color','red')
    xline(sigmatrue,'--','LineWidth',linesize)
    xline(sigmastart,'LineWidth',linesize,'Color','blue')
 %   xlim([0.2e-5 2e-5])
    title("\sigma")
end

%% True parameters
param_T_obs       = 7.8e-9;  % Reverberation time
param_G0_obs      = db2pow(-83.9);    % linear gain 
param_lambda_obs  = 10e9;       % arrival rate (1/s)    
param_sigma_N_obs       = sqrt(0.28e-9);   % Noise std

%% Prior ranges
 T_a = 10e-8; 
 T_b = 20e-9;  
 G0_a = db2pow(-90);    % Power gain (not in dB)
 G0_b = db2pow(-70);     % Power gain (not in dB)
 lambda_a = 1e8;
 lambda_b = 20e9;
 sigmaN_a = sqrt(2e-10); 
 sigmaN_b = sqrt(2e-9);
 
%% G0
figure
tiledlayout(4,N);

N = 10;
for i = 1:N
    nexttile
    [f,xi] = ksdensity(pow2db(params_G0(i,:)));
    plot(xi,f,'black')
    xlim([pow2db(G0_a) pow2db(G0_b)]) 
    if i == 1
        title('Iteration: 1')
    else
        title(i)
    end
    hold on
    % MAP estimate line
    [MAPest_val, index] = max(f);
    MAPline_xvals = ones(100,1)*xi(index);
    MAPline_yvals = linspace(0,MAPest_val,100);
    plot(MAPline_xvals,MAPline_yvals,'r')
    
    % True value line
    TRUEline_xvals = ones(100,1)*pow2db(param_G0_obs);
    TRUEline_yvals = linspace(0,MAPest_val,100);
    plot(TRUEline_xvals,TRUEline_yvals,'--g')
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    if i == 5
        xlabel("Initial gain, G_0")
    end
end
legend("Approx. posterior", "MAP estimate", "True value")

%% lambda
for i = 1:N
    nexttile
    [f,xi] = ksdensity(params_lambda(i,:));
    plot(xi,f,'black')
    xlim([lambda_a lambda_b])
    
    hold on
    % MAP estimate line
    [MAPest_val, index] = max(f);
    MAPline_xvals = ones(100,1)*xi(index);
    MAPline_yvals = linspace(0,MAPest_val,100);
    plot(MAPline_xvals,MAPline_yvals,'r')
    
    % True value line
    TRUEline_xvals = ones(100,1)*param_lambda_obs;
    TRUEline_yvals = linspace(0,MAPest_val,100);
    plot(TRUEline_xvals,TRUEline_yvals,'--g')
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    if i == 5
        xlabel("Arrival rate, \lambda")
    end
end

%% T
for i = 1:N
    nexttile
    [f,xi] = ksdensity(params_T(i,:));
    plot(xi,f,'black')
    xlim([T_a T_b])
    
        hold on
    % MAP estimate line
    [MAPest_val, index] = max(f);
    MAPline_xvals = ones(100,1)*xi(index);
    MAPline_yvals = linspace(0,MAPest_val,100);
    plot(MAPline_xvals,MAPline_yvals,'r')
    
    % True value line
    TRUEline_xvals = ones(100,1)*param_T_obs;
    TRUEline_yvals = linspace(0,MAPest_val,100);
    plot(TRUEline_xvals,TRUEline_yvals,'--g')
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    if i == 5
        xlabel("Reverberation time, T")
    end
end

%% sigma_N
for i = 1:N
    nexttile
    [f,xi] = ksdensity(params_sigma_N(i,:));
    plot(xi,f,'black')
    xlim([sigmaN_a sigmaN_b])
    
        hold on
    % MAP estimate line
    [MAPest_val, index] = max(f);
    MAPline_xvals = ones(100,1)*xi(index);
    MAPline_yvals = linspace(0,MAPest_val,100);
    plot(MAPline_xvals,MAPline_yvals,'r')
    
    % True value line
    TRUEline_xvals = ones(100,1)*param_sigma_N_obs;
    TRUEline_yvals = linspace(0,MAPest_val,100);
    plot(TRUEline_xvals,TRUEline_yvals,'--g')
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    if i == 5
        xlabel("Noise variance, \sigma_N^2")
    end
end
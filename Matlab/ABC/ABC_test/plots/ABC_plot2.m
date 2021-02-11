%% Plots
thetas = [params_T; params_G0; params_lambda; params_sigma_N];

tt = tiledlayout(num_iter,4,'TileSpacing','Compact','Padding','Compact');
linesize = 2;
truecolor = '#32CD32';

for i = 1:num_iter
    [f_T, xi_T] = ksdensity(params_T(i,:));
    [f_G0, xi_G0] = ksdensity(pow2db(params_G0(i,:)));
    [f_lambda, xi_lambda] = ksdensity(params_lambda(i,:));
    [f_sigma_N, xi_sigma_N] = ksdensity(params_sigma_N(i,:));
    
    nexttile
    area(xi_T*1e9,f_T,'FaceColor','#bbbbbb','LineStyle','None');
    MMSE_T = mean(params_T(i,:));
    xline(MMSE_T*1e9,'r','LineWidth',linesize)
    xline(theta_true(1)*1e9,'--','Color',truecolor,'LineWidth',linesize)
    subtitle('T \times 10^{-9}')
    xlim([prior(1,1)*1e9 prior(1,2)*1e9])
    %set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    if i == 1
        legend( "Approx. posterior",'MMSEE',"True value")
    end
    ylabel(['Iteration: ' num2str(i)])
    nexttile
    area(xi_G0,f_G0,'FaceColor','#bbbbbb','LineStyle','None')
    xline(pow2db(theta_true(2)),'--','Color',truecolor,'LineWidth',linesize)
    MMSE_G0 = mean(params_G0(i,:));
    xline(pow2db(MMSE_G0),'r','LineWidth',linesize)
    subtitle('G_0 [dB]')
    xlim([pow2db(prior(2,1)) pow2db(prior(2,2))])
    %set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    nexttile
    area(xi_lambda*1e-9,f_lambda,'FaceColor','#bbbbbb','LineStyle','None')
    xline(theta_true(3)*1e-9,'--','Color',truecolor,'LineWidth',linesize)
    MMSE_lambda = mean(params_lambda(i,:));
    xline(MMSE_lambda*1e-9,'r','LineWidth',linesize)
    subtitle('\lambda \times 10^{9}')
    xlim([prior(3,1)*1e-9 prior(3,2)*1e-9])
    %set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    nexttile
    area(xi_sigma_N.^2*1e9,f_sigma_N,'FaceColor','#bbbbbb','LineStyle','None')
    xline(theta_true(4)^2*1e9,'--','Color',truecolor,'LineWidth',linesize)
    MMSE_sigma_N = mean(params_sigma_N(i,:));
    xline(MMSE_sigma_N^2*1e9,'r','LineWidth',linesize)
    subtitle('\sigma_N^2 \times 10^{-9}')
    xlim([prior(4,1).^2*1e9 prior(4,2).^2*1e9])
    %set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
end
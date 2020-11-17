tiledlayout(2,2)
nexttile
h = area(xi_T,f_T,'FaceColor','#bbbbbb');
xline(theta_true(1),'--g')
xline(mean(thetas(1,:)),'r')
subtitle('T')
xlim([prior(1,1) prior(1,2)])
set(gca,'ytick',[])
set(gca,'yticklabel',[])


nexttile
area(xi_G0,f_G0,'FaceColor','#bbbbbb')
xline(theta_true(2),'--g')
xline(mean(thetas(2,:)),'r')
subtitle('G_0')
xlim([prior(2,1) prior(2,2)])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

nexttile
area(xi_lambda,f_lambda,'FaceColor','#bbbbbb')
xline(theta_true(3),'--g')
xline(mean(thetas(3,:)),'r')
subtitle('\lambda')
xlim([prior(3,1) prior(3,2)])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

nexttile
area(xi_sigma_N.^2,f_sigma_N,'FaceColor','#bbbbbb')
xline(theta_true(4)^2,'--g')
xline(mean(thetas(4,:))^2,'r')
subtitle('\sigma_N^2')
xlim([prior(4,1).^2 prior(4,2).^2])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
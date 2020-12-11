figure(1)
tiledlayout(3,1)
set(0,'defaulttextInterpreter','latex')
% for l = [1 4 9]
%    nexttile
%    plot(s_mean(1:i-1,l))
% end
i = 11;
   nexttile
   plot(s_sim(100:end,1))
   title('mean(m0)')
   set(gca,'ytick',[])
      nexttile
   plot(s_sim(100:end,4))
   title('Var(m1)')
   set(gca,'ytick',[])
      nexttile
   plot(s_sim(100:end,9))
   title('Cov(m0,m1)')
   set(gca,'ytick',[])
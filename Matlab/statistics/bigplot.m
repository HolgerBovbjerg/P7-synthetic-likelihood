%%
f2 = figure(2)
f2.Units = 'inches';
f2.OuterPosition = [1.25 1.25 10 14];

tmeanvar = tiledlayout(8,4)
%title(tmeanvar,"Mean and variance of the temporal moments",'Fontsize',24)
nexttile
plot(Tvary,mean(temporal0T,2),"o")
ylabel('Mean(m_0)','FontSize',12)

nexttile

plot(lambdavary,mean(temporal0lambda,2),"o")

nexttile
plot(sigma_noise_vary,mean(temporal0sigma_noise,2),"o")

nexttile
plot(G0vary,mean(temporal0G0,2),"o")


nexttile
plot(Tvary,mean(temporal1T,2),"o")
ylabel('Mean(m_1)','FontSize',12)

nexttile
plot(lambdavary,mean(temporal1lambda,2),"o")

nexttile
plot(sigma_noise_vary,mean(temporal1sigma_noise,2),"o")

nexttile
plot(G0vary,mean(temporal1G0,2),"o")

nexttile
plot(Tvary,mean(temporal2T,2),"o")
ylabel('Mean(m_2)','FontSize',12)

nexttile
plot(lambdavary,mean(temporal2lambda,2),"o")

nexttile
plot(sigma_noise_vary,mean(temporal2sigma_noise,2),"o")

nexttile
plot(G0vary,mean(temporal2G0,2),"o")

nexttile
plot(Tvary,mean(temporal3T,2),"o")
ylabel('Mean(m_3)','FontSize',12)

nexttile
plot(lambdavary,mean(temporal3lambda,2),"o")

nexttile
plot(sigma_noise_vary,mean(temporal3sigma_noise,2),"o")

nexttile
plot(G0vary,mean(temporal3G0,2),"o")




nexttile
plot(Tvary,var(temporal0T'),"o")
ylabel('Var(m_0)','FontSize',12)

nexttile
plot(lambdavary,var(temporal0lambda'),"o")

nexttile
plot(sigma_noise_vary,var(temporal0sigma_noise'),"o")

nexttile
plot(G0vary,var(temporal0G0'),"o")

nexttile
plot(Tvary,var(temporal1T'),"o")
ylabel('Var(m_1)','FontSize',12)

nexttile
plot(lambdavary,var(temporal1lambda'),"o")

nexttile
plot(sigma_noise_vary,var(temporal1sigma_noise'),"o")

nexttile
plot(G0vary,var(temporal1G0'),"o")

nexttile
plot(Tvary,var(temporal2T'),"o")
ylabel('Var(m_2)','FontSize',12)

nexttile
plot(lambdavary,var(temporal2lambda'),"o")

nexttile
plot(sigma_noise_vary,var(temporal2sigma_noise'),"o")

nexttile
plot(G0vary,var(temporal2G0'),"o")


nexttile
plot(Tvary,var(temporal3T'),"o")
xlabel("T",'FontSize',16)
ylabel('Var(m_3)','FontSize',12)

nexttile
plot(lambdavary,var(temporal3lambda'),"o")
xlabel("\lambda",'FontSize',16)

nexttile
plot(sigma_noise_vary,var(temporal3sigma_noise'),"o")
xlabel("\sigma_N",'FontSize',16)

nexttile
plot(G0vary,var(temporal3G0'),"o")
xlabel("G_0",'FontSize',16)

exportgraphics(gcf,'meanvarplot.pdf','ContentType','vector')

%%
f3 = figure(3)
f3.Units = 'inches';
f3.OuterPosition = [1.25 1.25 10 14];
lenmoment = 150;
tcorr = tiledlayout(6,4)
%title(tcorr,"Correlation between the temporal moments",'FontSize',24)
nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0T(qq,:),temporal1T(qq,:));
    pp(qq) = A(1,2);
end
plot(Tvary,pp,"o")
ylim([0 1])
ylabel('corr(m_0, m_1)','FontSize',12)

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0lambda(qq,:),temporal1lambda(qq,:));
    pp(qq) = A(1,2);
end
plot(lambdavary,pp,"o")
ylim([0 1])

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0sigma_noise(qq,:),temporal1sigma_noise(qq,:));
    pp(qq) = A(1,2);
end
plot(sigma_noise_vary,pp,"o")
ylim([0 1])


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0G0(qq,:),temporal1G0(qq,:));
    pp(qq) = A(1,2);
end
plot(G0vary,pp,"o")
ylim([0 1])

% ---------------------------------------------------------------------------------------------

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0T(qq,:),temporal2T(qq,:));
    pp(qq) = A(1,2);
end
plot(Tvary,pp,"o")
ylim([0 1])
ylabel('corr(m_0, m_2)','FontSize',12)


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0lambda(qq,:),temporal2lambda(qq,:));
    pp(qq) = A(1,2);
end
plot(lambdavary,pp,"o")
ylim([0 1])


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0sigma_noise(qq,:),temporal2sigma_noise(qq,:));
    pp(qq) = A(1,2);
end
plot(sigma_noise_vary,pp,"o")
ylim([0 1])

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0G0(qq,:),temporal2G0(qq,:));
    pp(qq) = A(1,2);
end
plot(G0vary,pp,"o")
ylim([0 1])


% ---------------------------------------------------------------------------------------------

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0T(qq,:),temporal3T(qq,:));
    pp(qq) = A(1,2);
end
plot(Tvary,pp,"o")
ylabel('corr(m_0, m_3)','FontSize',12)
ylim([0 1])


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0lambda(qq,:),temporal3lambda(qq,:));
    pp(qq) = A(1,2);
end
plot(lambdavary,pp,"o")
ylim([0 1])

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0sigma_noise(qq,:),temporal3sigma_noise(qq,:));
    pp(qq) = A(1,2);
end
plot(sigma_noise_vary,pp,"o")
ylim([0 1])


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal0G0(qq,:),temporal3G0(qq,:));
    pp(qq) = A(1,2);
end
plot(G0vary,pp,"o")
ylim([0 1])

% ---------------------------------------------------------------------------------------------

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal1T(qq,:),temporal2T(qq,:));
    pp(qq) = A(1,2);
end
plot(Tvary,pp,"o")
ylabel('corr(m_1, m_2)','FontSize',12)
ylim([0 1])


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal1lambda(qq,:),temporal2lambda(qq,:));
    pp(qq) = A(1,2);
end
plot(lambdavary,pp,"o")
ylim([0 1])


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal1sigma_noise(qq,:),temporal2sigma_noise(qq,:));
    pp(qq) = A(1,2);
end
plot(sigma_noise_vary,pp,"o")
ylim([0 1])

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal1G0(qq,:),temporal2G0(qq,:));
    pp(qq) = A(1,2);
end
plot(G0vary,pp,"o")
ylim([0 1])

% ---------------------------------------------------------------------------------------------

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal1T(qq,:),temporal3T(qq,:));
    pp(qq) = A(1,2);
end
plot(Tvary,pp,"o")
ylabel('corr(m_1, m_3)','FontSize',12)
ylim([0 1])


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal1lambda(qq,:),temporal3lambda(qq,:));
    pp(qq) = A(1,2);
end
plot(lambdavary,pp,"o")
ylim([0 1])

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal1sigma_noise(qq,:),temporal3sigma_noise(qq,:));
    pp(qq) = A(1,2);
end
plot(sigma_noise_vary,pp,"o")
ylim([0 1])

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal1G0(qq,:),temporal3G0(qq,:));
    pp(qq) = A(1,2);
end
plot(G0vary,pp,"o")
ylim([0 1])


% ---------------------------------------------------------------------------------------------

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal2T(qq,:),temporal3T(qq,:));
    pp(qq) = A(1,2);
end
plot(Tvary,pp,"o")
xlabel("T",'FontSize',16)
ylabel('corr(m_2, m_3)','FontSize',12)
ylim([0 1])


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal2lambda(qq,:),temporal3lambda(qq,:));
    pp(qq) = A(1,2);
end
plot(lambdavary,pp,"o")
xlabel("\lambda",'FontSize',16)
ylim([0 1])

nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal2sigma_noise(qq,:),temporal3sigma_noise(qq,:));
    pp(qq) = A(1,2);
end
plot(sigma_noise_vary,pp,"o")
xlabel("\sigma_N",'FontSize',16)
ylim([0 1])


nexttile
for qq = 1:lenmoment
    A = corrcoef(temporal2G0(qq,:),temporal3G0(qq,:));
    pp(qq) = A(1,2);
end
plot(G0vary,pp,"o")
xlabel("G_0",'FontSize',16)
ylim([0 1])


exportgraphics(f3,'correlationplot.pdf','ContentType','vector')
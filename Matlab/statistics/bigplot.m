%%

tiledlayout(8,4)
nexttile
plot(Tvary,mean(temporal0T,2),"o")

nexttile
plot(lambdavary,mean(temporal0lambda,2),"o")

nexttile
plot(sigma_noise_vary,mean(temporal0sigma_noise,2),"o")

nexttile
plot(G0vary,mean(temporal0G0,2),"o")


nexttile
plot(Tvary,mean(temporal1T,2),"o")

nexttile
plot(lambdavary,mean(temporal1lambda,2),"o")

nexttile
plot(sigma_noise_vary,mean(temporal1sigma_noise,2),"o")

nexttile
plot(G0vary,mean(temporal1G0,2),"o")

nexttile
plot(Tvary,mean(temporal2T,2),"o")

nexttile
plot(lambdavary,mean(temporal2lambda,2),"o")

nexttile
plot(sigma_noise_vary,mean(temporal2sigma_noise,2),"o")

nexttile
plot(G0vary,mean(temporal2G0,2),"o")

nexttile
plot(Tvary,mean(temporal3T,2),"o")

nexttile
plot(lambdavary,mean(temporal3lambda,2),"o")

nexttile
plot(sigma_noise_vary,mean(temporal3sigma_noise,2),"o")

nexttile
plot(G0vary,mean(temporal3G0,2),"o")




nexttile
plot(Tvary,var(temporal0T'),"o")

nexttile
plot(lambdavary,var(temporal0lambda'),"o")

nexttile
plot(sigma_noise_vary,var(temporal0sigma_noise'),"o")

nexttile
plot(G0vary,var(temporal0G0'),"o")

nexttile
plot(Tvary,var(temporal1T'),"o")

nexttile
plot(lambdavary,var(temporal1lambda'),"o")

nexttile
plot(sigma_noise,var(temporal1sigma_noise'),"o")

nexttile
plot(G0vary,var(temporal1G0'),"o")

nexttile
plot(Tvary,var(temporal2T'),"o")

nexttile
plot(lambdavary,var(temporal2lambda'),"o")

nexttile
plot(sigma_noise_vary,var(temporal2sigma_noise'),"o")

nexttile
plot(G0vary,var(temporal2G0'),"o")


nexttile
plot(Tvary,var(temporal3T'),"o")

nexttile
plot(lambdavary,var(temporal3lambda'),"o")

nexttile
plot(sigma_noise_vary,var(temporal3sigma_noise'),"o")

nexttile
plot(G0vary,var(temporal3G0'),"o")



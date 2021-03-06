%%
figure(6)
%covariance = covariancesmall;
normalT = makedist("Normal","mu",theta_true(1),"sigma",sqrt(covariance(1,1)))
normalG0 = makedist("Normal","mu",theta_true(2),"sigma",sqrt(covariance(2,2)))
normallambda = makedist("Normal","mu",theta_true(3),"sigma",sqrt(covariance(3,3)))
normalsigma = makedist("Normal","mu",theta_true(4),"sigma",sqrt(covariance(4,4)))


yT = linspace(prior(1,1),prior(1,2),500);
yG0 = linspace(prior(2,1),prior(2,2),500);
ylambda = linspace(prior(3,1),prior(3,2),500);
ysigma = linspace(prior(4,1),prior(4,2),500);

t = tiledlayout(4,1);
title(t,'Proposal distributions')

nexttile
plot(yT,pdf(normalT,yT))
xline(theta_true(1))
xlim([prior(1,1) prior(1,2)])

nexttile
plot(yG0,pdf(normalG0,yG0))
xline(theta_true(2))
xlim([prior(2,1) prior(2,2)])

nexttile
plot(ylambda,pdf(normallambda,ylambda))
xline(theta_true(3))
xlim([prior(3,1) prior(3,2)])

nexttile
plot(ysigma,pdf(normalsigma,ysigma))
xline(theta_true(4))
xlim([prior(4,1) prior(4,2)])

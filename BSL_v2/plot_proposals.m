%%
figure(4)

normalT = makedist("Normal","mu",theta_curr(1),"sigma",sqrt(covariance(1,1)))
normalG0 = makedist("Normal","mu",theta_curr(2),"sigma",sqrt(covariance(2,2)))
normallambda = makedist("Normal","mu",theta_curr(3),"sigma",sqrt(covariance(3,3)))
normalsigma = makedist("Normal","mu",theta_curr(4),"sigma",sqrt(covariance(4,4)))

yT = linspace(1e-8,1e-10,500);
yG0 = linspace(db2pow(-85),db2pow(-81),500);
ylambda = linspace(1e8,20e8,500);
ysigma = linspace(sqrt(0.1e-9),sqrt(0.4e-9),500);

tiledlayout(4,1)
nexttile
plot(yT,pdf(normalT,yT))

nexttile
plot(yG0,pdf(normalG0,yG0))

nexttile
plot(ylambda,pdf(normallambda,ylambda))

nexttile
plot(ysigma,pdf(normalsigma,ysigma))
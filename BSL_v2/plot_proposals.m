%%
figure(4)
% 
% normalT = makedist("Normal","mu",theta_curr(1),"sigma",sqrt(covariance(1,1)))
% normalG0 = makedist("Normal","mu",theta_curr(2),"sigma",sqrt(covariance(2,2)))
% normallambda = makedist("Normal","mu",theta_curr(3),"sigma",sqrt(covariance(3,3)))
% normalsigma = makedist("Normal","mu",theta_curr(4),"sigma",sqrt(covariance(4,4)))


normalT = makedist("Normal","mu",theta_curr(1),"sigma",sqrt(covariance(1,1)))
normalG0 = makedist("Normal","mu",theta_curr(2),"sigma",sqrt(covariance(2,2)))
normallambda = makedist("Normal","mu",theta_curr(3),"sigma",sqrt(covariance(3,3)))
normalsigma = makedist("Normal","mu",theta_curr(4),"sigma",sqrt(covariance(4,4)))


yT = linspace(0,2e-7,500);
yG0 = linspace(db2pow(-94),db2pow(-74),500);
ylambda = linspace(1e6,10e8,500);
ysigma = linspace(sqrt(0.28e-10),sqrt(0.28e-8),500);

tiledlayout(4,1)
nexttile
plot(yT,pdf(normalT,yT))

nexttile
plot(yG0,pdf(normalG0,yG0))

nexttile
plot(ylambda,pdf(normallambda,ylambda))

nexttile
plot(ysigma,pdf(normalsigma,ysigma))
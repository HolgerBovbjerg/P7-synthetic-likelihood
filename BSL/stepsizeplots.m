figure(5)

normalT = makedist("Normal","mu",7.8e-9,"sigma",1.958263692685272e-10)
normalG0 = makedist("Normal","mu",db2pow(-83.9),"sigma",9.115476637038641e-11)
normallambda = makedist("Normal","mu",10e8,"sigma",2.986250270870375e+07)
normalsigma = makedist("Normal","mu",sqrt(0.28e-9),"sigma",3.637616513134642e-07)

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
clear

tMin = 0;
tMax = 20;

lambda = 10;

numPoints = poissrnd(lambda);

t = tMax*rand(numPoints,1);
alpha = normrnd(0 ,1, numPoints,1);

scatter(t,alpha)
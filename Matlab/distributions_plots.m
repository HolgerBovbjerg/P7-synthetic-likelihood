%% Gaussian
clear 
x = (-5:0.001:5);
mu = 0;
sigma = 1;
fGauss = 1/(sigma*sqrt(2))*exp(-1/2*((x-mu)/sigma).^2);

figure
plot(x,fGauss, 'DisplayName',"\mu = " + mu + ", \sigma = " + sigma)
xlabel("x")
ylabel("f(x)")

hold on

mu = 1;
sigma = 1;
fGauss = 1/(sigma*sqrt(2))*exp(-1/2*((x-mu)/sigma).^2);
plot(x,fGauss, 'DisplayName',"\mu = " + mu + ", \sigma = " + sigma)



mu = 0;
sigma = 2;
fGauss = 1/(sigma*sqrt(2))*exp(-1/2*((x-mu)/sigma).^2);
plot(x,fGauss, 'DisplayName',"\mu = " + mu + ", \sigma = " + sigma)


lgd = legend;

lgd.Location = 'northwest';

exportgraphics(gcf,'gaussian_example.pdf','ContentType','vector')


%% Uniform
clear

x = -5:0.001:5;
b = 4;
a = -2;

fUnif = zeros(length(x),1);

for i = 1:length(x)
    if x(i) >= a && x(i) <= b 
        fUnif(i) = 1/(b-a);
    end
end

figure
plot(x,fUnif,'DisplayName', "a= "+a+", b= " + b)
hold on
lgd = legend;

b = 0;
a = -4;

fUnif = zeros(length(x),1);

for i = 1:length(x)
    if x(i) >= a && x(i) <= b 
        fUnif(i) = 1/(b-a);
    end
end

plot(x,fUnif,'DisplayName', "a= "+a+", b= " + b)
lgd = legend;
xlabel("x")
ylabel("f(x)")

b = 5;
a = 2;

fUnif = zeros(length(x),1);

for i = 1:length(x)
    if x(i) >= a && x(i) <= b 
        fUnif(i) = 1/(b-a);
    end
end

plot(x,fUnif,'DisplayName', "a= "+a+", b= " + b)

lgd = legend;

lgd.Location = 'northwest';

exportgraphics(gcf,'uniform_example.pdf','ContentType','vector')


%% Poisson
clear

k = 0:40;

lambda = 1;

fPoisson = lambda.^k*exp(-lambda)./factorial(k);

figure
plot(k,fPoisson,'-o', 'DisplayName', "\lambda= "+lambda)
hold on
xlabel("k")
ylabel("p(k)")

lambda = 3;

fPoisson = lambda.^k*exp(-lambda)./factorial(k);

plot(k,fPoisson,'-o', 'DisplayName', "\lambda= "+lambda)

lambda = 20;

fPoisson = lambda.^k*exp(-lambda)./factorial(k);

plot(k,fPoisson,'-o', 'DisplayName', "\lambda= "+lambda)


lgd = legend;

lgd.Location = 'northeast';

exportgraphics(gcf,'poisson_example.pdf','ContentType','vector')


%% Exponential
clear

x = 0:0.01:3;

lambda = 1;

fExponential = lambda.*exp(-lambda*x);
figure
plot(x,fExponential, 'DisplayName', "\lambda= "+lambda)
hold on
xlabel("x")
ylabel("f(x)")

lambda = 2;
fExponential = lambda.*exp(-lambda*x);

plot(x,fExponential, 'DisplayName', "\lambda= "+lambda)

lambda = 3;
fExponential = lambda.*exp(-lambda*x);

plot(x,fExponential, 'DisplayName', "\lambda= "+lambda)

lgd = legend;

lgd.Location = 'northeast';
xlim([0 2])
exportgraphics(gcf,'exponential_example.pdf','ContentType','vector')

%% Log-normal
clear 
x = (0:0.001:5);
mu = 0;
sigma = 1;
fLognorm = 1./(x*sigma*sqrt(2*pi)) .* exp(-1/2 * ((log(x)-mu)/(sigma)).^2);

figure
plot(x,fLognorm, 'DisplayName',"\mu = " + mu + ", \sigma = " + sigma)
xlabel("x")
ylabel("f(x)")

hold on

mu = 1;
sigma = 1;
fLognorm = 1./(x*sigma*sqrt(2*pi)) .* exp(-1/2 * ((log(x)-mu)/(sigma)).^2);
plot(x,fLognorm, 'DisplayName',"\mu = " + mu + ", \sigma = " + sigma)



mu = 1;
sigma = 0.5;
fLognorm = 1./(x*sigma*sqrt(2*pi)) .* exp(-1/2 * ((log(x)-mu)/(sigma)).^2);
plot(x,fLognorm, 'DisplayName',"\mu = " + mu + ", \sigma = " + sigma)

mu = 0.5;
sigma = 1;
fLognorm = 1./(x*sigma*sqrt(2*pi)) .* exp(-1/2 * ((log(x)-mu)/(sigma)).^2);
plot(x,fLognorm, 'DisplayName',"\mu = " + mu + ", \sigma = " + sigma)

lgd = legend;

lgd.Location = 'northeast';

exportgraphics(gcf,'lognormal_example.pdf','ContentType','vector')



%% Rayleigh
clear 
x = (0:0.001:5);

sigma = 0.5;
fRayleigh = x/(sigma^2).*exp(-1/2*((x/sigma).^2));

figure
plot(x,fRayleigh, 'DisplayName',"\sigma = " + sigma)
xlabel("x")
ylabel("f(x)")

hold on

sigma = 1;
fRayleigh = x/(sigma^2).*exp(-1/2*((x/sigma).^2));

plot(x,fRayleigh, 'DisplayName',"\sigma = " + sigma)

sigma = 2;
fRayleigh = x/(sigma^2).*exp(-1/2*((x/sigma).^2));

plot(x,fRayleigh, 'DisplayName',"\sigma = " + sigma)


lgd = legend;

lgd.Location = 'northeast';

exportgraphics(gcf,'rayleigh_example.pdf','ContentType','vector')
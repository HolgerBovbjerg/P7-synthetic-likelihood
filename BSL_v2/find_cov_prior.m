function [covariance thetacurr] = find_cov_prior()
Nl = 20000;

a = 9e-9; b = 2e-9;
T_prior = a + (b-a).*rand(Nl,1);
T = a + (b-a)*rand;

a = 5e8; b = 15e8;
lambda_prior = a + (b-a).*rand(Nl,1);
lambda = a + (b-a)*rand;


a = sqrt(0.1e-9); b = sqrt(0.5e-9);
sigma_prior = a + (b-a).*rand(Nl,1);
sigma = a + (b-a).*rand;


a = db2pow(-85); b = db2pow(-82);
G0prior = a + (b-a).*rand(Nl,1);
G0 = a + (b-a).*rand;

thetacurr = [T G0 lambda sigma]';
covariance =eye(4) .* [var(T_prior) var(G0prior) var(lambda_prior) var(sigma_prior)];


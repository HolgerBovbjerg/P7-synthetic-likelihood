function [covariance thetacurr] = find_cov_prior(param_G0_obs)
Nl = 20000;

T_min = 1e-9;
T_max = 15e-9;
T_prior = T_min + (T_max-T_min).*rand(Nl,1);
T = T_min + (T_max-T_min)*rand;

lambda_min = 1e8;
lambda_max = 20e9;
lambda_prior = lambda_min + (lambda_max-lambda_min).*rand(Nl,1);
lambda = lambda_min + (lambda_max-lambda_min)*rand;


sigmaN_min = sqrt(0.28e-10);
sigmaN_max = sqrt(0.28e-8);
sigma_prior = sigmaN_min + (sigmaN_max-sigmaN_min).*rand(Nl,1);
sigma = sigmaN_min + (sigmaN_max-sigmaN_min).*rand;

G0_min = db2pow(pow2db(param_G0_obs) - 10);    % Power gain (not in dB)
G0_max = db2pow(pow2db(param_G0_obs) + 10);     % Power gain (not in dB)
G0prior = G0_min + (G0_max-G0_min).*rand(Nl,1);
G0 = G0_min + (G0_max-G0_min).*rand;

thetacurr = [T G0 lambda sigma]';
covariance =eye(4) .* [var(T_prior) var(G0prior) var(lambda_prior) var(sigma_prior)];


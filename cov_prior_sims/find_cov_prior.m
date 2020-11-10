function [covariance, thetacurr] = find_cov_prior(prior)
Nl = 20000;

T_prior = prior(1,1) + (prior(1,2) - prior(1,1)).*rand(Nl,1);

G0_prior = prior(2,1) + (prior(2,2) - prior(2,1)).*rand(Nl,1);

lambda_prior = prior(3,1) + (prior(3,2) - prior(3,1)).*rand(Nl,1);

sigmaN_prior = prior(4,1) + (prior(4,2) - prior(4,1)).*rand(Nl,1);

thetacurr = [T_prior(1) G0_prior(1) lambda_prior(1) sigmaN_prior(1)]';
covariance = cov([T_prior G0_prior lambda_prior sigmaN_prior]);

% Use this for independent parameters
% covariance =eye(4) .* [var(T_prior) var(G0prior) var(lambda_prior) var(sigmaN_prior)];


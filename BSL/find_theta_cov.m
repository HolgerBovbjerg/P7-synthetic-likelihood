function covmatrix = find_theta_cov()
a = 9e-9; b = 2e-9,Nl = 2000;
Tvary = a + (b-a).*rand(Nl,1);

a = 5e8; b = 15e8;
lambdavary = a + (b-a).*rand(Nl,1);


a = sqrt(0.1e-9); b = sqrt(0.5e-9);
sigma_noise_vary = a + (b-a).*rand(Nl,1);


a = db2pow(-85); b = db2pow(-82);
G0vary = a + (b-a).*rand(Nl,1);

A = [sort(Tvary) sort(G0vary) sort(lambdavary) sort(sigma_noise_vary) ];
covmatrix = cov(A);

end
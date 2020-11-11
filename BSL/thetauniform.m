function theta = thetauniform()
a = 9e-9; b = 2e-9,Nl = 1;
T = a + (b-a).*rand(Nl,1);

a = 5e8; b = 15e8,Nl = 1;
lambda = a + (b-a).*rand(Nl,1)


a = sqrt(0.1e-9); b = sqrt(0.5e-9),Nl = 1;
sigma_noise = a + (b-a).*rand(Nl,1);


a = db2pow(-85); b = db2pow(-82),Nl = 1;
G0 = a + (b-a).*rand(Nl,1);
theta = [T lambda sigma_noise G0];
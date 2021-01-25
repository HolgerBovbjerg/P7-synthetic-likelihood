prior_T = [1e-9 100e-9]; 

prior_G0 =[db2pow(-94) db2pow(-74)];

prior_lambda = [1/200e-9 4e9];

prior_sigma = [sqrt(db2pow(-132)*800) sqrt(db2pow(-112)*800)];

prior = [prior_T; prior_G0; prior_lambda; prior_sigma];
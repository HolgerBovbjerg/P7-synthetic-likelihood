T_true = 7.8e-9;
G0_true = db2pow(-83.9); % Reverberation gain converted from dB to power - 4.0738e-9
lambda_true = 1e9; % Corresponding to 200 rays pr. 200ns
sigma_N_true = sqrt(0.28e-9); % noise standard deviation

theta_true = [T_true G0_true lambda_true sigma_N_true];
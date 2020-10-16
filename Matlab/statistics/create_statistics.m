% This function returns a matrix with N summary statistics vectors based
% upon N realizations from the Turin model. 
% Each realization is generated from parameters drawn uniformly between the
% limits given to the function. 
function [S T G0] = create_statistics(N, T_min, T_max, G0_min, G0_max, lambda_min, lambda_max, sigma_N_min, sigma_N_max)

    T = T_min + (T_max-T_min).*rand(N,1);
    G0 = G0_min + (G0_max-G0_min).*rand(N,1);
    lambda = lambda_min + (lambda_max-lambda_min).*rand(N,1);
    sigma_N = sigma_N_min + (sigma_N_max-sigma_N_min).*rand(N,1);
    
    B = 4e9;
    Ns = 801;
    
    S = zeros(N,4);
    for i = 1:N
    
        [Pv, t] = sim_turin_matrix_gpu(1, B, Ns, T(i), G0(i), lambda(i), sigma_N(i));

        t = gpuArray(t);
        S(i,1) = trapz(t,(t.^0.*Pv));
        S(i,2) = trapz(t,(t.^1.*Pv));
        S(i,3) = trapz(t,(t.^2.*Pv));
        S(i,4) = trapz(t,(t.^3.*Pv));

    end
end

% This function returns a matrix with N summary statistics vectors based
% upon N realizations from the Turin model. 
% Each realization is generated from parameters drawn uniformly between the
% limits given to the function. 
function S = create_statistics(N, T, G0, lambda, sigma_N, B, Ns)
    
    S = zeros(N,4);
    for i = 1:N
        % Draw a realization from the Turin model based on input parameters
        [Pv, t] = sim_turin_matrix_gpu(1, B, Ns, T, G0, lambda, sigma_N);

%         t = gpuArray(t); % Time it and see if it is faster with. 

        % Calculate temporal moments with numerical integration
        S(i,1) = trapz(t,(Pv));         % 0th moment (t^0*Pv = Pv)
        S(i,2) = trapz(t,(t.^1.*Pv));   % 1st moment
        S(i,3) = trapz(t,(t.^2.*Pv));   % 2nd moment
        S(i,4) = trapz(t,(t.^3.*Pv));   % 3rd moment

    end
end

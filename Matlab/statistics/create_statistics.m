% This function returns a matrix with N summary statistics vectors based
% upon N realizations from the Turin model. 
% Each realization is generated from parameters drawn uniformly between the
% limits given to the function. 
function S = create_statistics(M,N, T, G0, lambda, sigma_N, B, Ns)
    S = zeros(M,8);
    for j = 1:M
        % Draw a realization from the Turin model based on input parameters
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, T, G0, lambda, sigma_N);
        
        % t = gpuArray(t); % Time it and see if it is faster with. 
        % Calculate temporal moments with numerical integration
        m0 = trapz(t,(Pv));         % 0th moment (t^0*Pv = Pv)
        m1 = trapz(t,(t.^1.*Pv));   % 1st moment
        m2 = trapz(t,(t.^2.*Pv));   % 2nd moment
        m3 = trapz(t,(t.^3.*Pv));   % 3rd moment
       
        % Calculate summary statistics from moments
        S(j,1) = log(mean(m0));            % Mean of 0th moment 
        S(j,2) = log(mean(m1));            % Mean of 1st moment
        S(j,3) = log(mean(m2));            % Mean of 2nd moment
        S(j,4) = log(mean(m3));            % Mean of 3rd moment
        S(j,5) = log(var(m0));             % Variance of 0th moment
        S(j,6) = log(var(m1));             % Variance of 1st moment
        S(j,7) = log(var(m2));             % Variance of 2nd moment
        S(j,8) = log(var(m3));             % Variance of 3rd moment
    end
end

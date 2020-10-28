% This function returns a matrix with N summary statistics vectors based
% upon N realizations from the Turin model. 
% Each realization is generated from parameters drawn uniformly between the
% limits given to the function. 
function S = create_statistics(N, T, G0, lambda, sigma_N, B, Ns)
    
    S = zeros(1,4);
    
    % Draw a realization from the Turin model based on input parameters
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, T, G0, lambda, sigma_N);

    % t = gpuArray(t); % Time it and see if it is faster with. 
    % Calculate temporal moments with numerical integration
    m0 = trapz(t,(Pv));         % 0th moment (t^0*Pv = Pv)
    m1 = trapz(t,(t.^1.*Pv));   % 1st moment
%     m2 = trapz(t,(t.^2.*Pv));   % 2nd moment


    % Calculate summary statistics from moments
    S(1) = log(mean(m0));            % Mean of 0th moment 
    S(2) = log(mean(m1));            % Mean of 1st moment
%     S(3) = mean(m2);            % Mean of 2nd moment
    S(3) = log(var(m0));             % Variance of 0th moment
    S(4) = log(var(m1));             % Variance of 1st moment
%     S(6) = var(m2);             % Variance of 2nd moment


end

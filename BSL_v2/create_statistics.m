% This function returns a matrix with N summary statistics vectors based
% upon N realizations from the Turin model. 
% Each realization is generated from parameters drawn uniformly between the
% limits given to the function. 
function S = create_statistics(Pv, t)
    
    S = zeros(1,8);
    
    % Draw a realization from the Turin model based on input parameters
%     [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, T, G0, lambda, sigma_N);

    % t = gpuArray(t); % Time it and see if it is faster with. 

    % Calculate temporal moments with numerical integration
    m0_n2 = trapz(t,(Pv));           % 0th moment (t^0*Pv = Pv)
    m1_n2 = trapz(t,(t.^1.*Pv));     % 1st moment
    m0_n4 = trapz(t,(Pv).^2);         % 0th moment (t^0*Pv = Pv)
    m1_n4 = trapz(t,(t.^1.*Pv.^2));   % 1st moment    

%     % Calculate summary statistics from moments
%     S(1) = mean(log(m0));            % Mean of 0th moment 
%     S(2) = mean(log(m1));            % Mean of 1st moment
%     S(3) = mean(log(m2));            % Mean of 2nd moment
%     S(4) = mean(log(m3));            % Mean of 3rd moment
%     S(5) = var(log(m0));             % Variance of 0th moment
%     S(6) = var(log(m1));             % Variance of 1st moment
%     S(7) = var(log(m2));             % Variance of 2nd moment
%     S(8) = var(log(m3));             % Variance of 3rd moment
%     
        % Calculate summary statistics from moments
    S(1) = log(mean(m0_n2));            % Mean of 0th moment 
    S(2) = log(mean(m1_n2));            % Mean of 1st moment
    S(3) = log(mean(m0_n4));            % Mean of 2nd moment
    S(4) = log(mean(m1_n4));            % Mean of 3rd moment
    S(5) = log(var(m0_n2));             % Variance of 0th moment
    S(6) = log(var(m1_n2));             % Variance of 1st moment
    S(7) = log(var(m0_n4));             % Variance of 2nd moment
    S(8) = log(var(m1_n4));             % Variance of 3rd moment

end

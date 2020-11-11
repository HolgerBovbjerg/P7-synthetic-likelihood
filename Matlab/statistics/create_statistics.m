function S = create_statistics(Pv, t)
    S = zeros(1,4);
    % Calculate temporal moments with numerical integration
    m0 = trapz(t,(Pv));         % 0th moment (t^0*Pv = Pv)
    m1 = trapz(t,(t.^1.*Pv));   % 1st moment

    % Calculate summary statistics from moments
    S(1) = log(mean(m0));            % Ln of Mean of 0th moment 
    S(2) = log(mean(m1));            % Ln of Mean of 1st moment
    S(3) = log(var(m0));             % Ln of Variance of 0th moment
    S(4) = log(var(m1));             % Ln of Variance of 1st moment
end

function S = create_statistics_extended(Pv, t)
    S = zeros(1,18);
    % Calculate temporal moments with numerical integration
    m0 = trapz(t,Pv);         % 0th second order moment (t^0*Pv = Pv)
    m1 = trapz(t,t.^1.*Pv);   % 1st second order moment
    m2 = trapz(t,t.^2.*Pv);         % 2nd second order moment 
    m3 = trapz(t,t.^3.*Pv);   % 3rd second order moment
%     m0kurt = trapz(t,Pv);         % 0th fourth order moment (t^0*Pv = Pv)
%     m1kurt = trapz(t,t.^1.*(Pv.^2));   % 1st fourth order moment
%     m2kurt = trapz(t,t.^2.*(Pv.^2));         % 2nd fourth order moment 
%     m3kurt = trapz(t,t.^3.*(Pv.^2));   % 3rd fourth order moment

    % Calculate summary statistics from moments
    S(:) = [log(mean(m0));...            % Ln of Mean of 0th moment 
    log(mean(m1));...            % Ln of Mean of 1st moment
    log(mean(m2));...            % Ln of Mean of 2nd moment
    log(mean(m3));...            % Ln of Mean of 3rd moment
%     log(mean(m0kurt));...            % Ln of Mean of 0th moment
%     log(mean(m1kurt));...            % Ln of Mean of 1st moment
%     log(mean(m2kurt));...            % Ln of Mean of 2nd moment
%     log(mean(m3kurt));...            % Ln of Mean of 3rd moment
    log(var(m0));...             % Ln of Variance of 0th moment
    log(var(m1));...             % Ln of Variance of 1st moment
    log(var(m2));...             % Ln of Variance of 0th moment
    log(var(m3))];             % Ln of Variance of 1st moment
%     log(var(m0kurt));...             % Ln of Variance of 0th moment
%     log(var(m1kurt))];...             % Ln of Variance of 1st moment
%     log(var(m2kurt));...             % Ln of Variance of 0th moment
%     log(var(m3kurt))];             % Ln of Variance of 1st moment
end

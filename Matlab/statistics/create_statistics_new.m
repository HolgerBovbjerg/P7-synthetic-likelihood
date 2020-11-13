function S = create_statistics_new(Pv, t)
%     Pv = pow2db(Pv);
    % Calculate temporal moments with numerical integration
    m0 = trapz(t,Pv)';         % 0th second order moment (t^0*Pv = Pv)
    m1 = trapz(t,t.^1.*Pv)';   % 1st second order moment
    m2 = trapz(t,t.^2.*Pv)';         % 2nd second order moment 
%     m0kurt = trapz(t,Pv.^2)';         % 0th fourth order moment (t^0*Pv = Pv)
%     m1kurt = trapz(t,t.^1.*(Pv.^2))';   % 1st fourth order moment
%     m2kurt = trapz(t,t.^2.*(Pv.^2))';         % 2nd fourth order moment 

    D2 = [m0 m1 m2];
    D4 = [m0kurt m1kurt m2kurt];
    meanD2 = mean(D2);
%     meanD4 = mean(D4);
    CovD2 = cov(D2);
%     CovD4 = cov(D4);

    S = [meanD2(1) meanD2(2) meanD2(3)...
        CovD2(1,1) CovD2(2,2) CovD2(3,3)...
        CovD2(1,2) CovD2(1,3) CovD2(2,3)];%...
%   S = [meanD4(1) meanD4(2) meanD4(3)...
%         CovD4(1,1) CovD4(2,2) CovD4(3,3)...
%         CovD4(1,2) CovD4(1,3) CovD4(2,3)];
    S = log(S);
end

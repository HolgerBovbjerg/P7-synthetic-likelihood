function s_obs = observed_statistics
load('C:\Users\claus\Downloads\HfullNLOS5.mat')
Bw = 4e9;
Ns = 801;
deltaf = Bw/(Ns-1); % Frequency seperation
    tmax = 1/deltaf; % Maximum delay, found from the bandwidth via frequency seperation
    t = (0:Ns-1)'./(deltaf*Ns); % Generate timestamps, in seconds
    
P_h = zeros(50*50,801);
k = 1;
for l = 1:50
    for i = 1:50
        dat = squeeze(HfullNLOS5(i,l,:));  % Data = 801 entry vector with complex values. Corresponding to one antenna?
        P_h(k,:) = abs(ifft(dat,[],1)).^2; 
        k = k+1;
    end
end

s_obs = create_statistics(P_h', t);

end
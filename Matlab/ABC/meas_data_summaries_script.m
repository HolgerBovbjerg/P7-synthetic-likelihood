S_obs_meas_data = zeros(9,50);
delay = 26;
t = (0:800-delay)'*200e-9/(801-delay);

for j = 1:50
    Pv = zeros(50,801-delay);
    for i = 1:50
        Pv(i,:) = abs(ifft(squeeze(HfullNLOS5(j,i,delay+1:end)),[],2)).^2;
    end
    S_obs_meas_data(:,j) = create_statistics(Pv', t);
end


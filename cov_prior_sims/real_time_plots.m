function void = real_time_plots(thetas,j,accept,k)
T_true = 7.8e-9;
G0_true = db2pow(-83.9); % Reverberation gain converted from dB to power - 4.0738e-9
lambda_true = 10e8;
sigma_N_true = sqrt(0.28e-9); % noise variance

ttt = tiledlayout(4,1)
title(ttt,['Number of steps: ',num2str(j),'  Acceptance rate: ',num2str(accept/j)])
nexttile
plot(thetas(1,1:j),'o')
yline(T_true)
ylabel("Magnitude")
title("T")

nexttile
plot(pow2db(thetas(2,1:j)),'o')
yline(pow2db(G0_true))
ylabel("dB")
title("G0")

nexttile
plot(thetas(3,1:j),'o')
yline(lambda_true)
ylabel("Magnitude")
title("\lambda")

nexttile
plot(thetas(4,1:j),'o')
yline(sigma_N_true)
ylabel("Magnitude")
title("\sigma")

pause(.00000001);
if j ~=k
    delete(ttt)
end
void = 0;
end
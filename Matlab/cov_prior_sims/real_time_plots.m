function void = real_time_plots(theta_true, thetas, j, accept, k, prior)

ttt = tiledlayout(4,1);
title(ttt,['Number of steps: ',num2str(j),'  Acceptance rate: ',num2str(accept/j)])
nexttile
plot(thetas(1,1:j),'o')
ylim([prior(1,1) prior(1,2)])
yline(theta_true(1))
ylabel("Magnitude")
title("T")

nexttile
plot(pow2db(thetas(2,1:j)),'o')
ylim([pow2db(prior(2,1)) pow2db(prior(2,2))])
yline(pow2db(theta_true(2)))
ylabel("dB")
title("G0")

nexttile
plot(thetas(3,1:j),'o')
ylim([prior(3,1) prior(3,2)])
yline(theta_true(3))
ylabel("Magnitude")
title("\lambda")

nexttile
plot(thetas(4,1:j),'o')
ylim([prior(4,1) prior(4,2)])
yline(theta_true(4))
ylabel("Magnitude")
title("\sigma")

pause(.00000001);
if j ~=k-1
    delete(ttt)
end
void = 0;
end
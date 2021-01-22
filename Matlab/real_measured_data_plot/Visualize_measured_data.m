%% Ramoni code
Deltat = 1/diff([58 62]*1e9);
Taxis = (0:801-1)*Deltat;
% load ('HfullNLOS5.mat');

for ii = 1:50
    for jj = 1:50
        x(ii,jj,:) = ifft(squeeze(HfullNLOS5(ii,jj,:)));
    end
end

for ii=1:2
    for jj = 1:2
        if ii ==1 && jj ==1 
           h_vv = reshape(x(2*(ii-1)+1:2:end,2*(jj-1)+1:2:end,:),625,801);
        elseif ii == 1
           h_vh = reshape(x(2*(ii-1)+1:2:end,2*(jj-1):2:end,:),625,801);
        elseif jj == 1
            h_hv = reshape(x(2*(ii-1):2:end,2*(jj-1)+1:2:end,:),625,801);
        else
            h_hh = reshape(x(2*(ii-1):2:end,2*(jj-1):2:end,:),625,801);
        end
    end
end
Py_measured = abs(h_vv).^2;
Py_measured_ave = mean(Py_measured);


%% simulate data
load("Theta_true_values.mat");
Ns = 801;
B = 4e9;
deltaf = B/(Ns-1);

[Py_sim, ~] = sim_turin_matrix_gpu_w_delay(5000, B, Ns, theta_true, 6e-9);

%% create statistics
t = (0:Ns-1)'./(deltaf*Ns);

s_obs_measured = create_statistics(Py_measured',t);
s_obs_simulated = create_statistics(Py_sim,t);

%% plot data


figure
plot(t,pow2db( mean(abs(h_vv).^2) ))
hold on
plot(t,pow2db( mean(abs(h_vh).^2) ))
plot(t,pow2db( mean(abs(h_hv).^2) ))
plot(t,pow2db( mean(abs(h_hh).^2) ))

fig = figure;
plot(t*1e9,pow2db( mean(abs(h_vv).^2)) ,'LineWidth',2)
xtickformat('%d ns')
ylabel('[dB]')
box off

figure
plot(t,pow2db(Py_measured_ave),t,pow2db(mean(Py_sim,2)))


%% 
exportgraphics(fig,'observed_measurements_vertical_vertical_average_PDP.pdf','ContentType','vector')




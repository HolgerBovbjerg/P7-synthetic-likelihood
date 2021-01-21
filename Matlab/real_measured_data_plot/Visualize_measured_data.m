%%
clf
Deltat = 1/diff([58 62]*1e9);
Taxis = (0:801-1)*Deltat;
% load ('HfullNLOS5.mat');
for ii = 1:50
    for jj = 1:50
        x(ii,jj,:) =ifft(squeeze(HfullNLOS5(ii,jj,:)));
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


ift_hh = ifft(h_hh);
ift_vv = ifft(h_vv);
ift_hv = ifft(h_hv);
ift_vh = ifft(h_vh);




%% 
load('Theta_true_values.mat')

Theta_delay = [theta_true(1)/1.45 theta_true(2)/100e4 10e6 theta_true(4)/30];
Theta =       [theta_true(1)/1.45 theta_true(2)/10e4 10e6 theta_true(4)/30];
theta_test = theta_true;
theta_test(3) = 0.5e9; 
N = 300; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz

% Turin with delay
[P_y, t] = sim_turin_matrix_gpu_w_delay(N, B, Ns, Theta_delay, 1e-8);

[P_y2, t] = sim_turin_matrix_gpu(N, B, Ns, Theta);

[P_y3, t] = sim_turin_matrix_gpu(N, B, Ns, theta_test);
theta_test(3) = 2e9;
[P_y4, t] = sim_turin_matrix_gpu(N, B, Ns, theta_true);
[P_y5, t] = sim_turin_matrix_gpu(N, B, Ns, theta_test);
%%

plot(t, pow2db(mean(P_y3,2)), 'DisplayName','\lambda = 0.5e9')
hold on
plot(t, pow2db(mean(P_y4,2)), 'DisplayName','\lambda = 1e9')
plot(t, pow2db(mean(P_y5,2)), 'DisplayName','\lambda = 2e9')
legend

%%
t = (0:Ns-1)'./B;
figure(1) 
plot(t,pow2db(  abs(mean(ift_vv).^2) ) )

hold on 
plot(t, pow2db(mean(P_y,2)))
% plot(t, pow2db(mean(P_y2,2)))

% plot(t,pow2db(  abs(mean(ift_hh).^2) ) )
% % plot(t,pow2db(  abs(mean(ift_vh).^2) ) ) % Looks funny...
% plot(t,pow2db(  abs(mean(ift_hv).^2) ) )

%% 
exportgraphics(gcf,'observed_measurements_vertical_vertical_average_PDP.pdf','ContentType','vector')




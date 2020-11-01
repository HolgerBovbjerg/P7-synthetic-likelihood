%% G0
figure
tiledlayout(5,2)
for i = 1:10
    nexttile
    histogram(params_G0(i,:))
    xlim([db2pow(-95) db2pow(-78)]) 
end

%% lambda
figure
tiledlayout(5,2)
for i = 1:10
    nexttile
    histogram(params_lambda(i,:))
    xlim([1e9 1e11])
end

%% T
figure
tiledlayout(5,2)
for i = 1:10
    nexttile
    histogram(params_T(i,:))
    xlim([1e-10 1e-8])
end

%% sigma_N
figure
tiledlayout(5,2)
for i = 1:10
    nexttile
    histogram(params_sigma_N(i,:))
    xlim([sqrt(0.28e-10) sqrt(0.28e-8)])
end
clear 

%% 
data = HfullNLOS5;

Ns = length(data);
B = 4e9;

Yk = squeeze(data(3,1,:));

t = (0:Ns-1)'./B;
plot(t,pow2db(abs(ifft(Yk)).^2))
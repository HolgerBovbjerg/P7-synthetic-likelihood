
% Extract data from real-world measurement experiment:
% h_hh contains data from horizontal polarization on both receiving and transmitting antenna array 

% The resulting h_hh 625x801 matrix should be used to produce summary statistics 

Deltat = 1/diff([58 62]*1e9);
Taxis = (0:801-1)*Deltat;
%load ('HfullNLOS5.mat');
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
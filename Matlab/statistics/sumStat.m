
% Matlab script for finding summary statistics 

function [temporalMoment, temporalMoment_mean] = sumStat(t, data_set, nth_moment)
%% INPUT: 
%   data_set    = input vector n*m (row*column vector) with datapoints.
%   t           = input vector n*1 (row*column vector) with time stamps. 
%   nth_moment  = (integer) input the nth temporal moment we want out. 
%% OUTPUT:
%   temporalMoment        = vector with temporal moment.
%   temporalMoment_mean   = scalar with mean of all temporal moments 

len = length(data_set);
temporalMoment_0 = zeros(1,length(data_set));
data_set = data_set.'; % Transpose to rows

    for i = 1:len
        temporalMoment(i) = trapz(t,(t.^nth_moment.* data_set(:,i))); 
    end
    
temporalMoment       = transpose(temporalMoment);
temporalMoment_mean  = mean(temporalMoment);


%% OUTPUT: 


end



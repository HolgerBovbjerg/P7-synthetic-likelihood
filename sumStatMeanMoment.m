
% Matlab script for finding summary statistics based on moments  

function S = sumStatMeanMoment(t, data_set)
%% INPUT: 
%   data_set    = input vector n*m (row*column vector) with datapoints.
%   t           = input vector n*1 (row*column vector) with time stamps. 
%   nth_moment  = (integer) input the nth temporal moment we want out. 

%% OUTPUT:
%   S           = summary statistic vector with mean and variance of the 3 different moments. 

%% STEP 1: Find the temporal moments of the different turin realisations. 


sizeOne = size(data_set,2); % Returns the size of the column dimension 
                            % this is the number of Turin realizations 
 
% Prealocating vector for output summary statistics and vetors to hold temporal moments                          
S = zeros(6,1);                          
temporalMoment_0 = zeros(sizeOne,1); 
temporalMoment_1 = zeros(sizeOne,1); 
temporalMoment_2 = zeros(sizeOne,1);

% For the input data_set we have a matrix (row * column)
% Each column is one turin realization based on a specific set of parameters 
% Each row contains the value of the specific turin realization with a time stamp t.     

for i = 1:sizeOne  
    temporalMoment_0(i,1) = trapz(t,(t.^0.* data_set(:,i)));  % P_Y_obs(:,i) = the nth column 
    temporalMoment_1(i,1) = trapz(t,(t.^1.* data_set(:,i)));  % P_Y_obs(:,i) = the nth column 
    temporalMoment_2(i,1) = trapz(t,(t.^2.* data_set(:,i)));  % P_Y_obs(:,i) = the nth column 
     
end

% Create the summary statistic vector
S(1)  =  mean(temporalMoment_0);
S(2)  =  mean(temporalMoment_1);
S(3)  =  mean(temporalMoment_2);
S(4)  =  var(temporalMoment_0);
S(5)  =  var(temporalMoment_1);
S(6)  =  var(temporalMoment_2);

end



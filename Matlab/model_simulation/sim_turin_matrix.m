


function [P_Y, t] = sim_turin_matrix(N, Bw, Ns, T, G0, lambda, sigma_N)
    %% INPUTS:
    %  N       = Number of simulation iterations (The number of complete turin realizations)
    %  Bw      = Bandwidth 
    %  Ns      = Number of simulation iterations (The number of Hk realizations used) 
    %  T       = Reverberation time 
    %  G0      = Reverberation gain (power at delay zero)
    %  lambda  = Arrival rate of multipath components
    %  sigma_N = Noise variance
   
    %% OUTPUTS:
    % P_Y       = Power delay spectrum value - a mean of multiple touring simulations  
    % t         = Timestamp vector
 
    %% 
    deltaf = Bw/(Ns-1); % Frequency seperation
    tmax = 1/deltaf; % Maximum delay, found from the bandwidth via frequency seperation
    t = (0:Ns-1)'./(deltaf*Ns); % Generate timestamps, in seconds
    
    %% Simulate model
    % We then simulate the transfer function, H_k, of the channel, using the
    % Turin model. 
    Hk = zeros(Ns,N); % buffer for generated channel response data                                             
    % We run the simulation N times, creating new data sets for each
    % realization. 
    for n = 1:N
        lmax = 10000; % poissrnd(tmax*lambda);   % Number of multipath components, created from the Poisson distribution.
        tau = rand(lmax,1)*tmax;                 % time-delays, drawn uniformly between 0 and the maximum delay.  
        % For every multipath component a complex gain is generated, based on a
        % sigma generated from a corresponding delay time value. 
        % The complex number is generated in cartesian form by drawing the real
        % and the imaginary part seperately from a normal distribution. 
        sigma_alpha = sqrt(G0*exp(-(tau*(1/T)) ) / lambda);% Calculate variance using eq 13 and Ph = Lambda*sigma^2 from Ayush paper.
        % The complex valued alpha is created by combining the real and
        % imaginary parts.
        alpha = sigma_alpha * 1/sqrt(2) .* (randn(lmax,1) +  1j*randn(lmax,1));
        % For every frequency index, k, the contribution from every multipath
        % component is added to the transfer function. 
        k = (1:Ns);
        Hk(:,n) = (exp(-1j*2*pi*deltaf*k.*tau).' * alpha); % (Hk(:,n) - the nth column of matrix Hk)
    end
    
    % Generate noise vector complex normal distribution
    noise = sigma_N^2 * (1/sqrt(2))* (randn(Ns,N) + 1j*randn(Ns,N));
    % Add noise to transfer function
    Yk = Hk + noise;
    % Power delay profile:
    % Inverse furier transform of each column of Hk and then the absolute value
    % of this, elementwise squared. P_y -> A definition from Ayush paper.
    % P_y is now the Power delay profile in the time domain
    P_y = abs(ifft(Yk,[],1)).^2;
    % Averaging the power delay profile over the N realizations
    % P_y_mean = column vector where each entry is the mean of each row in
    % vector P_y
    P_y_mean = mean(P_y,2); 
    % We use P_Y = E_s * P_h + noise (Noise is already included in simulation)
    % E_s = energy of signal 
    P_Y = Bw*P_y_mean;
    
end


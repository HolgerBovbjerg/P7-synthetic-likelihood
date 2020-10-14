function [P_y, t] = sim_turin(N, B, Ns, T, G0, lambda, sigma_N)
    deltaf = B/(Ns-1); % Frequency seperation: 5 MHz
    tmax = 1/deltaf; % Maximum delay, found from the bandwidth via frequency seperation
    t = (0:Ns-1)'./(deltaf*Ns); % Generate timestamps, in seconds
    %% Simulate model
    % We then simulate the transfer function, H_k, of the channel, using the
    % Turin model. 

    Hk = zeros(Ns,N); % buffer for generated channel response data
    P_h = zeros(Ns,N); % Buffer for Power of impulse response
    ldist = makedist('poisson',tmax*lambda); % distribution for number 
                                                 % of multipaths to include is 
                                                 % a function of lambda and
                                                 % maximum delay
    % We run the simulation N times, creating new data sets for each
    % realization. 
    for n = 1:N
        lmax = random(ldist,1,1);   % Number of multipath components, created from the Poisson distribution.
        tau = rand(lmax,1)*tmax;    % time-delays, drawn uniformly between 0 and the maximum delay.  
        tau = sort(tau);
        % For every multipath component a complex gain is generated, based on a
        % sigma generated from a corresponding delay time value. 
        % The complex number is generated in cartesian form by drawing the real
        % and the imaginary part seperately from a normal distribution. 
        for i = 1:lmax 
                sigma_alpha = sqrt(G0*exp(-(tau(i)/T) ) / lambda);% Calculate variance using eq 13 and Ph = Lambda*sigma^2 from Ayush paper.
                % The complex valued alpha is created by combining the real and
                % imaginary parts. 
                alpha = 1/sqrt(2)*sigma_alpha*(randn +  1j*randn);  
                % For every frequency index, k, the contribution from multipath
                % component, i, is added to the transfer function. 
                for k = 1:Ns 
                    Hk(k,n) = Hk(k,n) + alpha*exp(-1j*2*pi*deltaf*k*tau(i));
                end
        end
        %% Generate noise vector
        noise = sigma_N^2 * (1/sqrt(2))* (randn(Ns,1) + 1j*randn(Ns,1));
        %% Add noise to transfer function with complex normal dist
        Hk(:,n) = Hk(:,n) + noise;
        % Power delay profile
        P_h(:,n) = abs(ifft(Hk(:,n))).^2;
    end
    %% Averaging over the N realizations
    P_h_mean = mean(P_h,2); % acg. power delay profile

    %% Simulating the power spectrum, P_h, from the transfer function
    % We use P_Y = E_s * P_h + noise (Noise is already included in simulation)
    P_y = B*P_h_mean; %+ sigma_N^2/Ns;
end
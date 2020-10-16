function L = synth_likelihood(s,Nr)
    mu = sum(s./N);
    S = s-mu;
    Sigma = S*S'./(Nr-1);
    L = -1/2*((s-mu)'\Sigma)*S - 1/2*log10(Sigma);
end
